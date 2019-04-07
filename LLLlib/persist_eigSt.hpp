
/*!
  \file persist_eigSt.hpp
*/

#ifndef PERSIST_EIGST_HPP
#define PERSIST_EIGST_HPP

#include <LLLlib/eigSt.hpp>
#include <iostream>
#include <cassert>


/*! \class persist_eigSt : public eigSt
  \brief Supplies the ability of persisting the state's data by inserting into a stream
  or extracting it from a stream
*/
class persist_eigSt : public eigSt {

  bool m_bBinary;

public:

  ////////////////////////////////
 // constructors / destructors //
////////////////////////////////

  /*!
    \fn persist_eigSt ();
    \brief constucts empty states
 
    persist_eigSt ();

    \fn persist_eigSt (Basis *basis_to_use, bool bComplex = false)
    \brief state to store WF
    \param basis_to_use: basis to be used for this state
    \param bComplex: wether coefficients are stored as complex
  */
  persist_eigSt (Basis *basis_to_use, bool bComplex = false);
  

  /*!
    \fn persist_eigSt (const persist_eigSt &src)
    \brief copy-constructor
  */
  persist_eigSt (const persist_eigSt &src);

  /*!
    \fn void set_binary(bool yes)
    \brief set binary-mode. 
    \param yes : true -> stream insertor will insert persist_eigSt as binary data into stream; 
    false -> stream insertor will insert persist_eigSt as ascii data into stream
  */
  void set_binary(bool yes);

  /*!
    \fn bool is_binary (void) const
    \brief retrieves, wether persist_eigSt will be saved as binary (true) or as ascii (false)
  */
  bool is_binary (void) const { return m_bBinary; };

  /*!
    \fn ~persist_eigSt ();
    \brief destructor
  */
  virtual ~persist_eigSt ();


  ////////////////////
 // public members //
////////////////////

  /*
    \fn void reinitialize(Basis *pBasis_to_use)
    \brief initialization (i.e. after array-allocation...)
  
  void reinitialize(Basis *pBasis_to_use);
  */

  //////////////////////
 // friend Operators //
//////////////////////

  /*!
    \fn operator friend std::ostream & operator<< (std::ostream &out, const persist_eigSt &eig)
    \brief stream insertor
    inserts data of eigSt into stream depending on mode in m_bbinary    
  */
  friend std::ostream & operator<< (std::ostream &out, const persist_eigSt &eig);

  /*!
    \fn friend std::istream & operator>> (std::istream &in, persist_eigSt &eig)
    \brief stream extractor
    extracts data of eigSt from stream depending on mode in m_bbinary
  */
  friend std::istream & operator>> (std::istream &in, persist_eigSt &eig);

};


  //////////////////////
 // freind operators //
//////////////////////

  // stream inserter
  std::ostream & operator<< (std::ostream &out, const persist_eigSt &eigst);

  // stream extracter
  std::istream & operator>> (std::istream &in, persist_eigSt &eigst);

#endif
