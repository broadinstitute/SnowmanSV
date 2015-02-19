//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ReadTable - A 0-indexed table of reads
//
#include <iostream>
#include <algorithm>
#include "ReadTable.h"
#include "SeqReader.h"
//#include "contigs.h"

// Read the sequences from a file
ReadTable::ReadTable(std::string filename, uint32_t reader_flags)
{
    idx = 0;
    m_pIndex = NULL; // not built by default
    SeqReader reader(filename, reader_flags);
    SeqRecord sr;
    while(reader.get(sr))
    {
        addRead(sr.toSeqItem());
    }
}

// JEREMIAH
/*ReadTable::ReadTable(std::string* seq, std::string * id, int length)
{
    m_pIndex = NULL; // not built by default
    for (int i = 0; i < length; i++) {
      SeqItem r;
      r.id = id[i];
      std::cout << id[i] << "\n";
      r.seq = seq[i];
      m_table.push_back(r);
    }
    }*/

// JEREMIAH
ReadTable::ReadTable(SeqRecordVector srv)
{
    idx = 0;
    m_pIndex = NULL; // not built by default
    for(std::vector<SeqRecord>::size_type i = 0; i != srv.size(); i++) {
      m_table.push_back(srv[i].toSeqItem());
    }
}

// JEREMIAH
ReadTable::ReadTable(const ContigVector &contigs)
{
    idx = 0;
    m_pIndex = NULL; // not built by default
    for (ContigVector::const_iterator it = contigs.begin(); it != contigs.end(); it++) {
      SeqItem si;
      si.seq = it->getSeq();
      si.id = it->getID();
      m_table.push_back(si);
    }

}

// JEREMIAH
ReadTable::ReadTable(const BamAlignmentVector &bav)
{
    idx = 0;
    m_pIndex = NULL; // not built by default
    for(BamAlignmentVector::const_iterator i = bav.begin(); i != bav.end(); i++) { 
      SeqItem si;
      if (!i->GetTag("SR", si.id)) {
	cout << "Expecting SR tag (unique read name)" << endl;
	exit(EXIT_FAILURE);
      }
      //si.id = i->Name;
      //si.seq = i->QueryBases;
      std::string seqr;

      if (!i->GetTag("TS", seqr))
	seqr = i->QueryBases;
      si.seq = seqr;
      assert(seqr.length());
      m_table.push_back(si);
    }
}

// JEREMIAH
ReadTable::ReadTable(const BamAlignmentUPVector &bav)
{
    idx = 0;
    m_pIndex = NULL; // not built by default
    for (auto &i : bav) {
      //for(BamAlignmentUPVector::const_iterator i = bav.begin(); i != bav.end(); i++) { 
      SeqItem si;
      if(!i->GetTag("SR", si.id)) {
	cerr << "ERROR: expecting SR tag for reads" << endl;
	exit(EXIT_FAILURE);
      } 

      std::string seqr;

      // if quality trimmed bases not found, go with original
      if (!i->GetTag("TS", seqr))
	seqr = i->QueryBases;
      si.seq = seqr;
      assert(seqr.length());
      m_table.push_back(si);
    }
}

// 
ReadTable::~ReadTable()
{
    if(m_pIndex != NULL)
        delete m_pIndex;
}

// Populate this read table with the reverse reads from pRT
void ReadTable::initializeReverse(const ReadTable* pRT)
{
    size_t numReads = pRT->getCount();
    m_table.reserve(numReads);
    for(size_t i = 0; i < numReads; ++i)
    {
        SeqItem read = pRT->getRead(i);
        read.seq.reverse();
        addRead(read);
    }
}

//
void ReadTable::reverseAll()
{
    for(size_t i = 0; i < getCount(); ++i)
        m_table[i].seq.reverse();
}

//
void ReadTable::addRead(const SeqItem& r)
{
    m_table.push_back(r);
}

//
size_t ReadTable::getReadLength(size_t idx) const
{
    return m_table[idx].seq.length();
}

//
std::string ReadTable::getReadID(size_t idx) const
{
  return m_table[idx].id;
}


//
const SeqItem& ReadTable::getRead(size_t idx) const
{
    assert(idx < m_table.size());
    return m_table[idx];
}

// Jeremiah
bool ReadTable::getRead(SeqItem &si)
{
  if (idx >= m_table.size())
    return false;
  si = m_table[idx];
  idx++;
  return true;
}

//JEREMIAH 
void ReadTable::setZero() 
{
  idx = 0;
} 


// indexReadsByID must be called before this function can be used
const SeqItem& ReadTable::getRead(const std::string& id) const
{
    assert(m_pIndex != NULL);
    if(m_pIndex == NULL)
    {
        std::cerr << "Error: read table is not indexed (did you forget to call ReadTable::buildIndex?)\n";
        assert(false);
    }

    ReadIndex::const_iterator i = m_pIndex->find(id);
    if(i == m_pIndex->end())
    {
        std::cerr << "Read with id " << id << " not found in table\n";
        assert(false);
    }

    return *i->second;
}

// build a read id -> *seqitem index
void ReadTable::indexReadsByID()
{
    m_pIndex = new ReadIndex;
    for(size_t i = 0; i < m_table.size(); ++i)
    {
        std::pair<ReadIndex::iterator, bool> result = 
            m_pIndex->insert(std::make_pair(m_table[i].id, &m_table[i]));
        if(!result.second)
        {
            std::cerr << "Error: Attempted to insert vertex into table with a duplicate id: " <<
                         m_table[i].id << "\n";
            std::cerr << "All reads must have a unique identifier\n";
            exit(1);
        }
    }
}

//
size_t ReadTable::getCount() const
{
    return m_table.size();
}

// 
size_t ReadTable::countSumLengths() const
{
    size_t sum = 0;
    for(size_t i = 0; i < m_table.size(); ++i)
    {
        sum += m_table[i].seq.length();
    }
    return sum;
}

//
void ReadTable::clear()
{
    m_table.clear();

    if(m_pIndex != NULL)
    {
        delete m_pIndex;
        m_pIndex = NULL;
    }
}

//
std::ostream& operator<<(std::ostream& out, const ReadTable& rt)
{
    size_t numReads = rt.getCount();
    for(size_t i = 0; i < numReads; ++i)
    {
        const SeqItem& read = rt.getRead(i);
        out << read.id << "\t" << read.seq.toString() << "\n";
    }
    return out;
}
