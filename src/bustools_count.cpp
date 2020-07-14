#include <iostream>
#include <fstream>
#include <cstring>

#include "Common.hpp"
#include "BUSData.h"

#include "bustools_count.h"


void bustools_count(Bustools_opt &opt) {
  BUSHeader h;
  size_t nr = 0;
  size_t N = 100000;
  uint32_t bclen = 0;
  BUSData* p = new BUSData[N];
  bool ignore_umi = false;

  // read and parse the equivelence class files

  std::unordered_map<std::vector<int32_t>, int32_t, SortedVectorHasher> ecmapinv;
  std::vector<std::vector<int32_t>> ecmap;

  std::unordered_map<std::string, int32_t> txnames;
  parseTranscripts(opt.count_txp, txnames);
  std::vector<int32_t> genemap(txnames.size(), -1);
  std::unordered_map<std::string, int32_t> genenames;
  parseGenes(opt.count_genes, txnames, genemap, genenames);
  parseECs(opt.count_ecs, h);
  ecmap = std::move(h.ecs);
  ecmapinv.reserve(ecmap.size());
  for (int32_t ec = 0; ec < ecmap.size(); ec++) {
    ecmapinv.insert({ecmap[ec], ec});
  }
  std::vector<std::vector<int32_t>> ec2genes;        
  create_ec2genes(ecmap, genemap, ec2genes);


  std::ofstream of;
  std::string mtx_ofn = opt.output + ".mtx";
  std::string barcodes_ofn = opt.output + ".barcodes.txt";
  std::string ec_ofn = opt.output + ".ec.txt";
  std::string gene_ofn = opt.output + ".genes.txt";
  of.open(mtx_ofn); 

  // write out the initial header
  of << "%%MatrixMarket matrix coordinate real general\n%\n";
  // number of genes
  auto mat_header_pos = of.tellp();
  std::string dummy_header(66, '\n');
  for (int i = 0; i < 33; i++) {
    dummy_header[2*i] = '%';
  }
  of.write(dummy_header.c_str(), dummy_header.size());


  size_t n_cols = 0;
  size_t n_rows = 0;
  size_t n_entries = 0;
  std::vector<BUSData> v;
  v.reserve(N);
  uint64_t current_bc = 0xFFFFFFFFFFFFFFFFULL;
  //temporary data
  std::vector<int32_t> ecs;
  std::vector<int32_t> glist;
  ecs.reserve(100);
  std::vector<int32_t> u;
  u.reserve(100);
  std::vector<int32_t> column_v;
  std::vector<std::pair<int32_t, double>> column_vp;
  if (!opt.count_collapse) {
    column_v.reserve(N); 
  } else {
    column_vp.reserve(N);
    glist.reserve(100);
  }
  //barcodes 
  std::vector<uint64_t> barcodes;
  int bad_count = 0;
  int compacted = 0;
  int rescued = 0;

  auto write_barcode_matrix = [&](const std::vector<BUSData> &v) {
    if(v.empty()) {
      return;
    }
    column_v.resize(0);
    n_rows+= 1;
    
    barcodes.push_back(v[0].barcode);
    double val = 0.0;
    size_t n = v.size();

    for (size_t i = 0; i < n; ) {
      size_t j = i+1;
      for (; j < n; j++) {
        if (v[i].UMI != v[j].UMI || ignore_umi) {
          break;
        }
      }

      // v[i..j-1] share the same UMI
      ecs.resize(0);
      for (size_t k = i; k < j; k++) {
        ecs.push_back(v[k].ec);
      }

      int32_t ec = v[i].ec;/*intersect_ecs(ecs, u, genemap, ecmap, ecmapinv, ec2genes);
      if (ec == -1) {
        ec = intersect_ecs_with_genes(ecs, genemap, ecmap, ecmapinv, ec2genes);              
        if (ec == -1) {
          bad_count += j-i;
        } else {
          bool filter = false;
          if (!opt.count_gene_multimapping) {
            filter = (ec2genes[ec].size() != 1);
          }
          if (!filter) {
            rescued += j-i;
            column_v.push_back(ec);
          }
        }

      } else {
        bool filter = false;
        if (!opt.count_gene_multimapping) {
          filter = (ec2genes[ec].size() != 1);
        }
        if (!filter) {
          compacted += j-i-1;
          column_v.push_back(ec);
        }
      }*/
      column_v.push_back(ec);
      i = j; // increment
    }
    std::sort(column_v.begin(), column_v.end());
    size_t m = column_v.size();
    for (size_t i = 0; i < m; ) {
      size_t j = i+1;
      for (; j < m; j++) {
        if (column_v[i] != column_v[j]) {
          break;
        }
      }
      double val = j-i;
      if (ignore_umi) {
        val = v[i].count;
      }
      of << n_rows << " " << (column_v[i]+1) << " " << val << "\n";
      n_entries++;
      
      i = j; // increment
    }
  };

  auto write_barcode_matrix_collapsed = [&](const std::vector<BUSData> &v) {
    if(v.empty()) {
      return;
    }
    column_vp.resize(0);
    n_rows+= 1;
    
    barcodes.push_back(v[0].barcode);
    double val = 0.0;
    size_t n = v.size();

    std::vector<std::vector<int32_t>> ambiguous_genes;

    for (size_t i = 0; i < n; ) {
      size_t j = i+1;
      for (; j < n; j++) {
        if (v[i].UMI != v[j].UMI || ignore_umi) {
          break;
        }
      }

      // v[i..j-1] share the same UMI
      ecs.resize(0);
      for (size_t k = i; k < j; k++) {
        ecs.push_back(v[k].ec);
      }

      intersect_genes_of_ecs(ecs,ec2genes, glist);
      int gn = glist.size();
      if (gn > 0) {
        if (opt.count_gene_multimapping) {
          for (auto x : glist) {
            column_vp.push_back({x, 1.0/gn});
          }
        } else {
          if (gn==1) {
            column_vp.push_back({glist[0],1.0});
          } else if (opt.count_em) {
            ambiguous_genes.push_back(std::move(glist));
          }
        }
      }
      i = j; // increment
    }
    std::sort(column_vp.begin(), column_vp.end());
    size_t m = column_vp.size();
    std::unordered_map<int32_t, double> col_map(m);
    std::vector<int32_t> cols;

    for (size_t i = 0; i < m; ) {
      size_t j = i+1;
      double val = column_vp[i].second;
      for (; j < m; j++) {
        if (column_vp[i].first != column_vp[j].first) {
          break;
        }
        val += column_vp[j].second;
      }
      col_map.insert({column_vp[i].first,val});
      cols.push_back(column_vp[i].first);

      n_entries++;
      
      i = j; // increment
    }

    if (opt.count_em) {
      //std::cerr << "Running EM algorithm" << std::endl;
      std::unordered_map<int32_t, double> c1,c2;
      // initialize with unique counts
      for (const auto &x : cols) {
        double val = 0;
        auto it = col_map.find(x);
        if (it != col_map.end()) {
          val = it->second;
        }
        c1.insert({x,val});
        c2.insert({x,0.0});
      }

      // c1 is from previous round, c2 is the new one
      int round = 0;
      double tol = 1e-3;
      double err = 1.0;
      if (!ambiguous_genes.empty()) {        
        while (round < 100 && err >= tol) {
          for (auto &x : col_map) {
            c2[x.first] = x.second;
          }

          for (const auto &glist : ambiguous_genes) {
            double s = 0.0;
            for (const auto &x : glist) {
              auto it = c1.find(x);
              if (it != c1.end()) {
                s += it->second;
              }
            }
            if (s > 0.0) {
              for (const auto &x : glist) {
                auto it = c1.find(x);
                if (it != c1.end()) {
                  c2[x] += it->second/s;
                }
              }
            }
          }

          err = 0.0;
          for (const auto &x : c1) {
            double a=x.second,b=c2[x.first];
            if (a>0 && b > 0) {
              err += std::abs((a-b)/(a+b));
            }
          }
          //std::cerr << "Round " << round << " error = " << err << std::endl;
          std::swap(c1,c2);
          round++;
        }
        std::swap(col_map,c1);
      }

    }



    for (const auto &x : cols) {
      double val = 0;
      auto it = col_map.find(x);
      if (it != col_map.end()) {
        val = it->second;
      }
      of << n_rows << " " << (x+1) << " " << val << "\n";
    }
  };

  for (const auto& infn : opt.files) { 
    std::streambuf *inbuf;
    std::ifstream inf;
    if (!opt.stream_in) {
      inf.open(infn.c_str(), std::ios::binary);
      inbuf = inf.rdbuf();
    } else {
      inbuf = std::cin.rdbuf();
    }
    std::istream in(inbuf); 

    parseHeader(in, h);
    bclen = h.bclen;
    if (h.umilen == 0) {
      ignore_umi = true;
    }
    
    int rc = 0;
    while (true) {
      in.read((char*)p, N*sizeof(BUSData));
      size_t rc = in.gcount() / sizeof(BUSData);
      nr += rc;
      if (rc == 0) {
        break;
      }

      
      for (size_t i = 0; i < rc; i++) {
        if (p[i].barcode != current_bc) {                 
          // output whatever is in v
          if (!v.empty()) {
            if (!opt.count_collapse) {
              write_barcode_matrix(v);
            } else {
              write_barcode_matrix_collapsed(v);
            }
          }
          v.clear();
          current_bc = p[i].barcode;
        }
        v.push_back(p[i]);

      }            
    }
    if (!v.empty()) {
      if (!opt.count_collapse) {
        write_barcode_matrix(v);
      } else {
        write_barcode_matrix_collapsed(v);
      }
    }

    if (!opt.stream_in) {
      inf.close();
    }
  }
  delete[] p; p = nullptr;

  if (!opt.count_collapse) {
    n_cols = ecmap.size();
  } else {
    n_cols = genenames.size();
  }

  of.close();
  
  std::stringstream ss;
  ss << n_rows << " " << n_cols << " " << n_entries << "\n";
  std::string header = ss.str();
  int hlen = header.size();
  assert(hlen < 66);
  of.open(mtx_ofn, std::ios::binary | std::ios::in | std::ios::out);
  of.seekp(mat_header_pos);
  of.write("%",1);
  of.write(std::string(66-hlen-2,' ').c_str(),66-hlen-2);
  of.write("\n",1);
  of.write(header.c_str(), hlen);
  of.close();

  // write updated ec file
  h.ecs = std::move(ecmap);
  if (!opt.count_collapse) {
    writeECs(ec_ofn, h);
  } else {
    writeGenes(gene_ofn, genenames);
  }
  // write barcode file
  std::ofstream bcof;
  bcof.open(barcodes_ofn);
  for (const auto &x : barcodes) {
    bcof << binaryToString(x, bclen) << "\n";
  }
  bcof.close();
  //std::cerr << "bad counts = " << bad_count <<", rescued  =" << rescued << ", compacted = " << compacted << std::endl;

  //std::cerr << "Read in " << nr << " BUS records" << std::endl;
}
