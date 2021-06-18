#include <iostream>
#include <fstream>
#include <thread>
#include <mutex>
#include <zlib.h>
#include <unordered_set>

#include "kseq.h"
#include "edlib.h"
#include "Common.hpp"
#include "BUSData.h"

#include "bustools_tempo.h"

KSEQ_INIT(gzFile, gzread);

class FastqSequenceReader {
public:
  FastqSequenceReader(const Bustools_opt& opt, int batch_barcode = 0, bool check_batch_barcode = false) : 
  state(false), readbatch_id(-1), current_file(0), numreads(0), batch_barcode(batch_barcode), check_batch_barcode(check_batch_barcode) {
    files.clear();
    for (int i = 0; i < opt.fastq.size(); ++i) {
      files.push_back(opt.fastq[i]);
    }
    nfiles = opt.nFastqs;
    reserveNfiles(nfiles);
  }
  FastqSequenceReader() : state(false), readbatch_id(-1), current_file(0), numreads(0), batch_barcode(0), check_batch_barcode(false) {};
  FastqSequenceReader(FastqSequenceReader &&o);
  ~FastqSequenceReader();
  
  bool empty();  
  void reset();
  void reserveNfiles(int n);
  size_t fetchSequences(char *buf, const int limit, std::vector<std::pair<const char*, int>>& seqs,
                      std::vector<std::pair<const char*, int>>& names,
                      std::vector<std::pair<const char*, int>>& quals,
                      int &readbatch_id,
                      bool full,
                      BUSData* p,
                      size_t bus_i_start,
                      size_t rc);
  
public:
  int nfiles = 1;
  std::vector<gzFile> fp;
  std::vector<int> l;
  std::vector<int> nl;
  std::vector<std::string> files;
  int current_file;
  std::vector<kseq_t*> seq;
  bool state; // is the file open
  int readbatch_id;
  int batch_barcode;
  bool check_batch_barcode;
  uint64_t numreads;
};

FastqSequenceReader::~FastqSequenceReader() {
  for (auto &f : fp) {
    if (f) {
      gzclose(f);
    }
  }
  for (auto &s : seq) {
    kseq_destroy(s);
  }
}

bool FastqSequenceReader::empty() {
  return (!state && current_file >= files.size());
}

void FastqSequenceReader::reset() {
  state = false;
  readbatch_id = -1;
  for (auto &f : fp) {
    if (f) {
      gzclose(f);
    }
    f = nullptr;
  }
  
  for (auto &ll : l) {
    ll = 0;
  }
  for (auto &nll : nl) {
    nll = 0;
  }
  
  current_file = 0;
  for (auto &s : seq) {
    kseq_destroy(s);
    s = nullptr;
  }
}

void FastqSequenceReader::reserveNfiles(int n) {
  fp.resize(nfiles);
  l.resize(nfiles, 0);
  nl.resize(nfiles, 0);
  seq.resize(nfiles, nullptr);
}

// returns true if there is more left to read from the files
size_t FastqSequenceReader::fetchSequences(char *buf, const int limit, std::vector<std::pair<const char *, int> > &seqs,
                                         std::vector<std::pair<const char *, int> > &names,
                                         std::vector<std::pair<const char *, int> > &quals,
                                         int& read_id,
                                         bool full, BUSData* p, size_t bus_i_start, size_t rc) {
  assert(bus_i_start < rc);
  std::string line;
  readbatch_id += 1; // increase the batch id
  read_id = readbatch_id; // copy now because we are inside a lock
  seqs.clear();
  if (full) {
    names.clear();
    quals.clear();
  }
  int bufpos = 0;
  int pad = nfiles;
  size_t bus_i = bus_i_start;
  while (true) {
    if (!state) { // should we open a file
      if (current_file >= files.size()) {
        // nothing left
        return 0;
      } else {
        // close the current files
        for (auto &f : fp) {
          if (f) {
            gzclose(f);
          }
        }
        // open the next one
        for (int i = 0; i < nfiles; i++) {
          fp[i] = gzopen(files[current_file+i].c_str(), "r");
          seq[i] = kseq_init(fp[i]);
          l[i] = kseq_read(seq[i]);
          
        }
        current_file+=nfiles;
        state = true; 
      }
    }
    // the file is open and we have read into seq1 and seq2
    bool all_l = true;
    int bufadd = nfiles;
    for (int i = 0; i < nfiles; i++) {
      all_l = all_l && l[i] >= 0;
      bufadd += l[i]; // includes seq
    }
    if (all_l) {      
      // fits into the buffer
      if (full) {
        for (int i = 0; i < nfiles; i++) {
          nl[i] = seq[i]->name.l;
          bufadd += l[i] + nl[i]; // includes name and qual
        }
        bufadd += 2*pad;
      }

      if (bufpos+bufadd< limit && bus_i < rc) {
        if (check_batch_barcode && p[bus_i].barcode != batch_barcode) {
          bus_i++;
          continue; // skip this BUS record
        } else if (p[bus_i].flags != numreads) {
          numreads++;
          continue; // skip this read
        } else if (p[bus_i].flags < numreads) {
          std::cerr << "Error: BUS file not sorted by flag; Encountered flag: " << p[bus_i].flags << 
            " while trying to read fastq record number: "<< numreads << std::endl;
          exit(1);
        }
        for (int i = 0; i < nfiles; i++) {
          char *pi = buf + bufpos;
          //std::cout << seq[i]->seq.s << std::endl; // DEBUG
          memcpy(pi, seq[i]->seq.s, l[i]+1);
          bufpos += l[i]+1;
          seqs.emplace_back(pi,l[i]);
          if (full) {
            pi = buf + bufpos;
            memcpy(pi, seq[i]->qual.s,l[i]+1);
            bufpos += l[i]+1;
            quals.emplace_back(pi,l[i]);
            pi = buf + bufpos;
            memcpy(pi, seq[i]->name.s, nl[i]+1);
            bufpos += nl[i]+1;
            names.emplace_back(pi, nl[i]);
          }
        }
        numreads++;
        bus_i++;
      } else {
        return rc - bus_i; // Return how much there's still left to read from the BUS records currently supplied
      }
      
      // read for the next one
      for (int i = 0; i < nfiles; i++) {
        l[i] = kseq_read(seq[i]);
      }        
    } else {
      state = false; // haven't opened file yet
    }
  }
}

FastqSequenceReader::FastqSequenceReader(FastqSequenceReader&& o) :
  nfiles(o.nfiles),
  batch_barcode(o.batch_barcode),
  check_batch_barcode(o.check_batch_barcode),
  fp(std::move(o.fp)),
  l(std::move(o.l)),
  nl(std::move(o.nl)),
  files(std::move(o.files)),
  current_file(o.current_file),
  numreads(o.numreads),
  seq(std::move(o.seq)) {
  
  o.fp.resize(nfiles);
  o.l.resize(nfiles, 0);
  o.nl.resize(nfiles, 0);
  o.seq.resize(nfiles, nullptr);
  o.state = false;
}

class TempoMasterProcessor {
public:
  TempoMasterProcessor(const Bustools_opt& opt, int nbatches) : 
  opt(opt), numreads(0), numreadsprocessed(0), N(100000), 
  nr(0), rc(0), nbatches(nbatches), finished_reading(false), in(NULL) {
    if (nbatches == 1) {
      FastqSequenceReader fSR(opt);
      FSRs.push_back(std::move(fSR));
    } else {
      for (int i = 0; i < nbatches; i++) {
        FastqSequenceReader fSR(opt, i, true); // Note: For batch mode, barcodes in BUS file must have values 0,1,2,...
        fSR.files.erase(fSR.files.begin(), fSR.files.begin()+opt.nFastqs*i);
        fSR.files.erase(fSR.files.begin()+opt.nFastqs, fSR.files.end());
        assert(fSR.files.size() == opt.nFastqs);
        FSRs.push_back(std::move(fSR));
      }
    }
    if (!opt.stream_in) {
      inf.open(opt.files[0].c_str(), std::ios::binary);
      inbuf = inf.rdbuf();
    } else {
      inbuf = std::cin.rdbuf();
    }
    in.rdbuf(inbuf);
    parseHeader(in, h);
    parseECs(opt.count_ecs, h);
    parseTranscripts(opt.count_txp, txnames);
    if (!parseFasta(opt.fasta, txnames, txseqs)) {
      exit(1);
    }
    p = new BUSData[N];
  }
  
  ~TempoMasterProcessor() {
    delete[] p;
  }
  
  BUSHeader h;
  Bustools_opt opt;
  int nbatches;
  std::unordered_set<int> curr_batches;
  std::unordered_set<int> curr_threads;
  uint64_t numreads;
  uint64_t numreadsprocessed;
  size_t N;
  size_t nr;
  size_t rc;
  std::mutex read_mutex;
  std::mutex write_mutex;
  bool finished_reading;
  std::streambuf *inbuf;
  std::ifstream inf;
  std::istream in;
  BUSData *p;
  std::unordered_map<std::string, int32_t> txnames;
  std::unordered_map<int32_t, std::pair <std::string, std::string> > txseqs;
  std::vector<std::mutex> fastq_mutexes;
  std::vector<FastqSequenceReader> FSRs;
  
  bool getBUSData(int i, int id) {
    std::lock_guard<std::mutex> lock(read_mutex);
    if (finished_reading) {
      return false;
    }
    curr_batches.emplace(i);
    curr_threads.emplace(id);
    if (curr_batches.size() >= nbatches && curr_threads.size() >= opt.threads) {
      in.read((char *) p, N * sizeof(BUSData));
      rc = in.gcount() / sizeof(BUSData);
      nr += rc;
      curr_batches.clear();
      curr_threads.clear();
      if (rc == 0) {
        finished_reading = true;
        return false;
      }
      return true;
    }
    return false;
  }
  
  
  
  bool doneReading() {
    return finished_reading;
  }
  
  void update(uint64_t n, uint64_t new_total_n) {
    std::lock_guard<std::mutex> lock(write_mutex);
    numreadsprocessed += n;
    numreads = new_total_n;
  }  
};

class TempoReadProcessor {
public:
  TempoReadProcessor(const Bustools_opt& opt, 
                    TempoMasterProcessor& mp,
                    int nfiles,
                    int nbatches,
                    int id) 
    : bufsize(1ULL<<23), 
      paired (opt.fastq_paired),
      opt (opt), 
      mp (mp),
      nfiles(nfiles),
      nbatches (nbatches),
      id (id) {
        buffer = new char[bufsize];
        seqs.reserve(bufsize/50);
    }
  
  TempoReadProcessor(TempoReadProcessor && o)
    : bufsize (o.bufsize),
      paired (o.paired),
      opt (o.opt),
      mp (o.mp),
      nfiles (o.nfiles),
      nbatches (o.nbatches),
      id (o.id),
      seqs (std::move(o.seqs)),
      names (std::move(o.names)),
      quals (std::move(o.quals)),
      bv(std::move(o.bv))
    { buffer = o.buffer; o.buffer = nullptr; o.bufsize = 0; }
      
  ~TempoReadProcessor() {
    if (buffer != nullptr) {
      delete[] buffer;
      buffer = nullptr;
    }
  }
  
  char *buffer;
  size_t bufsize;
  int nfiles = 1;
  int nbatches;
  int id;
  bool paired;
  Bustools_opt opt;
  TempoMasterProcessor& mp;
  std::vector<std::pair<const char*, int>> seqs;
  std::vector<std::pair<const char*, int>> names;
  std::vector<std::pair<const char*, int>> quals;
  std::vector<BUSData> bv;
  
  void operator()() {
    uint64_t read_counter = 0; // loop counter
    while (true) {
      int readbatch_id;
      int i = read_counter % nbatches;
      read_counter++;
      if (mp.getBUSData(i,id)) {
        size_t bus_i_start = 0;
        size_t remaining;
        uint64_t numreads;
        uint64_t n;
        do {
          {
          std::lock_guard<std::mutex> lock(mp.fastq_mutexes[i]);
          remaining = mp.FSRs[i].fetchSequences(buffer, bufsize, seqs, names, quals, readbatch_id, false, mp.p, bus_i_start, mp.rc);
          bus_i_start = mp.rc - remaining;
          numreads = mp.FSRs[i].numreads;
          n = seqs.size() / mp.FSRs[i].nfiles;
          }
          processBuffer();
          mp.update(n, numreads);
          clear();
        } while (remaining != 0);
      } else if (mp.doneReading()) {
        return;
      }
    }
  }
  
  void processBuffer() {
    int incf, jmax;
    incf = nfiles-1;
    jmax = nfiles;
    
    std::vector<const char*> s(jmax, nullptr);
    std::vector<int> l(jmax,0);
    
    for (int i = 0; i + incf < seqs.size(); i++) {
      for (int j = 0; j < jmax; j++) {
        s[j] = seqs[i+j].first;
        l[j] = seqs[i+j].second;      
      }
      i += incf;
      
      // find where the sequence is
      const char *seq = nullptr;
      const char *seq2 = nullptr;
      size_t seqlen = 0;
      size_t seqlen2 = 0;
      
      if (!paired) {
        seq = s[0]; // Read the 0th file
        seqlen = l[0];
      }
    }
  }
  
  void clear() {
    memset(buffer,0,bufsize);
    bv.clear();
  }
};

void bustools_tempo(const Bustools_opt& opt) {
  std::vector<std::thread> workers;
  assert(opt.fastq.size() % opt.nFastqs == 0);
  int nbatches = 1; // Normally, only need one sequence reader (except in the case of batch mode)
  if (!opt.batch_file.empty()) {
    nbatches = opt.fastq.size() / opt.nFastqs; // flag column starts at 0 for each run if BUS file generated in batch mode
  }
  TempoMasterProcessor mp(opt, nbatches);
  std::vector<std::mutex> mutexes(nbatches);
  mp.fastq_mutexes.swap(mutexes);
  
  for (int i = 0; i < opt.threads; i++) {
    workers.emplace_back(std::thread(TempoReadProcessor(opt,mp,opt.nFastqs,nbatches,i)));
  }
  
  // let the workers do their thing
  for (int i = 0; i < opt.threads; i++) {
    workers[i].join(); //wait for them to finish
  }
  
  std::cout << "Read in " << mp.nr << " BUS records" << std::endl;
  if (mp.numreadsprocessed < mp.nr) {
    std::cerr << "Warning: number of reads in FASTQs was less than number of reads in BUS file" << std::endl;
  }
  std::cout << "Read a total of " << mp.numreads << " reads from FASTQs, of which " 
  << mp.numreadsprocessed << " reads were further processed" << std::endl;
}



