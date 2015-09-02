#include "BamSplitter.h"

#include "htslib/khash.h"

void BamSplitter::splitBam() {

  assert(m_writers.size() == m_frac.size());

  std::vector<double> csum;
  double cs = 0;
  for (auto& i : m_frac) {
    cs += i;
    csum.push_back(cs);
  }
    
  SnowTools::BamRead r;
  
  bool rule_pass;

  std::cerr << "**Starting to split BAM " << __startMessage(); 

  size_t countr = 0;
  std::vector<size_t> all_counts(m_writers.size(), 0);
  while (GetNextRead(r, rule_pass))
    {

      // print a message
      ++countr;
      if (countr % 1000000 == 0) 
	std::cerr << "...at position " << r.Brief() << std::endl; 

      // put in one BAM or another
      uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(r.Qname().c_str()) ^ m_seed);
      double hash_val = (double)(k&0xffffff) / 0x1000000;

      for (size_t i = 0; i < m_writers.size(); ++i) {
	// decide whether to keep
	if (hash_val <= csum[i]/*sample_rate*/) {
	  r.RemoveAllTags();
	  m_writers[i].WriteAlignment(r);
	  ++all_counts[i];
	  break;
	}
      }
      
      
    }
  
}

void BamSplitter::fractionateBam(const std::string& outbam, SnowTools::Fractions& f) {

  std::cerr << "...setting up output fractionated BAM " << outbam << std::endl;
  
  SnowTools::BamWalker w;
  bam_hdr_t * h = bam_hdr_dup(br.get());
  
  w.SetWriteHeader(h);
  w.OpenWriteBam(outbam);

  f.m_frc.createTreeMap();

  SnowTools::BamRead r;
  
  bool rule_pass;
  
  std::cerr << "**Starting to fractionate BAM " << __startMessage(); 

  size_t countr = 0;
  std::vector<size_t> all_counts(m_writers.size(), 0);
  while (GetNextRead(r, rule_pass))
    {
      
      // print a message
      ++countr;
      if (countr % 1000000 == 0) 
	std::cerr << "...at position " << r.Brief() << std::endl; 
      
      SnowTools::GenomicRegion gr(r.ChrID(), r.Position(), r.Position(), '*');
      std::vector<int32_t> qid, sid;
      SnowTools::GRC grc(gr);
      grc.createTreeMap();
      f.m_frc.findOverlaps(grc, qid, sid, true);

      /*
      std::cerr << gr << " qid.size() " << qid.size() << " sid.size() " << sid.size() << std::endl;
      for (auto& i : qid) 
	std::cerr << "       Q: " << i << std::endl;
      for (auto& i : sid) 
	std::cerr << "       S: " << i << std::endl;
      */

      if (!qid.size())
	continue;
      
      double sample_rate = f.m_frc.at(qid[0]).frac;
      
      // put in one BAM or another
      uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(r.Qname().c_str()) ^ m_seed);
      double hash_val = (double)(k&0xffffff) / 0x1000000;
      
      //std::cerr << "sample rate " << sample_rate << " hash " << hash_val << std::endl;

      // decide whether to keep
      if (hash_val <= sample_rate) {
	//r.RemoveAllTags();
	w.WriteAlignment(r);
      }
    }
}

void BamSplitter::setWriters(const std::vector<std::string>& writers, const std::vector<double>& fracs) {

  assert(fracs.size() == writers.size());  

  for (auto& i : writers) {

    std::cerr << "...setting up BAM: " << i << std::endl;
    SnowTools::BamWalker w;
    bam_hdr_t * h = bam_hdr_dup(br.get()/*header()*/);
    w.SetWriteHeader(h);
    w.OpenWriteBam(i);
    m_writers.push_back(w);

  }

  m_frac = fracs;


}

std::string BamSplitter::__startMessage() const {

  std::stringstream ss;
  if (m_region.size() == 1)
    ss << " on region " << m_region[0] << std::endl;
  else if (m_region.size() == 0)
    ss << " on WHOLE GENOME" << std::endl;
  else
    ss << " on " << m_region.size() << " regions " << std::endl;
  return ss.str();

}