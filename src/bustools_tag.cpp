#include <iostream>
#include <fstream>
#include <algorithm>

#include "Common.hpp"
#include "BUSData.h"

#include "bustools_tag.h"

void bustools_tag(Bustools_opt &opt) {  
  uint32_t f = 0;
  uint32_t hamming_dist = 1;
  uint64_t tag = stringToBinary(opt.tag, f);
  uint32_t umilen = 0;
  uint32_t taglen = (opt.tag).size();
  
  BUSHeader h;
  bool outheader_written = false;
  
  std::streambuf *buf = nullptr;
  std::ofstream of;
  of.open(opt.output, std::ios::out | std::ios::binary);
  buf = of.rdbuf();
  std::ostream bus_out(buf);
  
  size_t nr = 0, nt = 0;
  size_t N = 100000;
  BUSData* p = new BUSData[N];
  BUSData bd;
  
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

    if (!outheader_written) {
      writeHeader(bus_out, h);
      outheader_written = true;
    }

    if (umilen == 0) {
      umilen = h.umilen;
    }
    
    if (umilen > taglen) {
      std::cerr << "Error: Tag length of " << taglen << " cannot be greater than UMI length of " << umilen << std::endl;
      exit(1);
    }

    int rc = 0;
    while (true) {
      in.read((char*)p, N*sizeof(BUSData));
      size_t rc = in.gcount() / sizeof(BUSData);
      if (rc == 0) {
        break;
      }
      nr +=rc;

      for (size_t i = 0; i < rc; i++) {
        bd = p[i];
        if (hamming(stringToBinary(tag, f), bd.UMI >> 2*(umi_len-taglen), taglen) <= hamming_dist) {
          bd.flags = 0;
          bus_out.write((char*) &bd, sizeof(bd));
          ++nt;
        } else {
          bd.flags = 1;
          bus_out.write((char*) &bd, sizeof(bd));    
        }
      }
    }
    if (!opt.stream_in) {
      inf.close();
    }
  }
  
  if (of.is_open()) {
    of.close();
  }
  
  delete[] p; p = nullptr;
  
  std::cerr << "Read in " << nr << " BUS records" << std::endl
    << "Tagged = " << nt << std::endl;
}
