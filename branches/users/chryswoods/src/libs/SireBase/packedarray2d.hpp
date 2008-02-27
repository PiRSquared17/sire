/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the developer's mailing list
  *  at http://siremol.org
  *
\*********************************************/

#ifndef SIREBASE_PACKEDARRAY2D_HPP
#define SIREBASE_PACKEDARRAY2D_HPP

#include "packedarray2d.h"
#include <boost/assert.hpp>

SIRE_BEGIN_HEADER

namespace SireBase
{

template<class T>
class PackedArray2D;

namespace detail
{

template<class T>
class SharedArray2DPtr;

template<class T>
class PackedArray2DData;

template<class T>
class PackedArray2D_ArrayData;

template<class T>
class PackedArray2D_Array;

/** This is a simple class that is used to manage the memory and 
    reference counting for the PackedArray2D */
template<class T>
class PackedArray2DMemory : public PackedArray2DMemoryBase
{
public:
    static char* create(quint32 narrays, quint32 nvalues);
    static char* detach(char *this_ptr, quint32 this_idx);
    static void destroy(PackedArray2DData<T> *arrays);
    
    static void incref(char *this_ptr, quint32 this_idx);
    static void decref(char *this_ptr, quint32 this_idx);

    static SharedArray2DPtr< PackedArray2DData<T> > shared_null;
};

/** This is the implicitly shared pointer class that 
    is used to hold any of the PackedArray2DMemory allocated objects
    
    @author Christopher Woods
*/
template<class T>
class SharedArray2DPtr
{

public:
    SharedArray2DPtr() : ptr(0)
    {}

    SharedArray2DPtr(const T *p)
    {
        ptr = const_cast<T*>(p);
    
        if (ptr)
            ptr->incref();
    }
    
    SharedArray2DPtr(const SharedArray2DPtr &other) : ptr(other.ptr)
    {
        if (ptr)
            ptr->incref();
    }
    
    ~SharedArray2DPtr()
    {
        if (ptr)
            ptr->decref();
    }
    
    SharedArray2DPtr<T>& operator=(const SharedArray2DPtr &other)
    {
        if (ptr != other.ptr)
        {
            T *new_ptr = other.ptr;
            
            //increment the other reference count
            if (new_ptr)
                new_ptr->incref();
                
            //decrement our reference count
            if (ptr)
                ptr->decref();
                
            //set the new pointer
            ptr = new_ptr;
        }
        
        return *this;
    }
    
    const T& operator*() const
    {
        return *ptr;
    }
    
    const T* operator->() const
    {
        return ptr;
    }
    
    T& operator*()
    {
        if (ptr)
            ptr = ptr->detach();
            
        return *ptr;
    }
    
    T* operator->()
    {
        if (ptr)
            ptr = ptr->detach();
            
        return ptr;
    }
    
    const T* data() const
    {
        return ptr;
    }
    
    const T* constData() const
    {
        return ptr;
    }
    
    T* data()
    {
        if (ptr)
            ptr = ptr->detach();
        
        return ptr;
    }
    
    /** Assign this pointer to point at 'weakptr' 
        but *without* changing the reference count.
        You ABSOLUTELY MUST ensure that you call 
        SharedArray2DPtr::weakRelease() before this 
        pointer is deleted or reassigned, so 
        as to not decrement the reference count incorrectly! */
    void weakAssign(T *weakptr)
    {
        if (ptr)
            ptr->decref();
            
        ptr = weakptr;
    }
    
    /** Release the pointer *without* decrementing the
        reference count. You should only call this
        function if the pointer was assigned using 
        the 'weakAssign()' function */
    void weakRelease()
    {
        ptr = 0;
    }
    
private:
    /** Actual pointer */
    T *ptr;
};

/** This class is used to hold all of the metadata about the 
    packed array of arrays
*/
template<class T>
class PackedArray2DData : public PackedArray2DDataBase
{

friend class PackedArray2DMemory<T>;

public:
    PackedArray2DData();
    PackedArray2DData(quint32 narrays, quint32 nvalues);
    
    PackedArray2DData(const PackedArray2DData &other);
    
    ~PackedArray2DData();

    void decref();
    
    PackedArray2DData* detach();
    
    const PackedArray2D_ArrayData<T>* nullArray() const;
    
    const PackedArray2D_ArrayData<T>* arrayDataData() const;
    const PackedArray2D_Array<T>* arrayData() const;

    const T* valueData() const;

    PackedArray2D_ArrayData<T>* arrayDataData();
    PackedArray2D_Array<T>* arrayData();

    T* valueData();
    
    void setNValuesInArray(quint32 i, quint32 nvalues);
    
    void close();
};

/** This class holds the metadata about an individual array
    in the packed collection of arrays */
template<class T>
class PackedArray2D_ArrayData : public PackedArray2D_ArrayDataBase
{

friend class PackedArray2D<T>;
friend class PackedArray2DData<T>;
friend class PackedArray2DMemory<T>;

public:
    PackedArray2D_ArrayData();
    PackedArray2D_ArrayData(quint32 this_idx);
    
    PackedArray2D_ArrayData(const PackedArray2D_ArrayData &other);
    
    ~PackedArray2D_ArrayData();
    
    void incref();
    void decref();
    
    PackedArray2D_ArrayData<T>* detach();
    
    PackedArray2DData<T>* extract() const;
    
    const T* valueData() const;
    T* valueData();
};

/** This is the implicitly shared array class that is packed
    together to make a PackedArray2D<T>. This class is mainly used
    via its alias - PackedArray2D<T>::Array. Note that the size
    of a PackedArray2D<T>::Array is fixed at construction. You
    cannot change the size as this array is packed contiguously
    in memory with other arrays.
    
    @author Christopher Woods
*/
template<class T>
class PackedArray2D_Array
{

friend class PackedArray2DMemory<T>;
friend class PackedArray2D<T>;

public:
    PackedArray2D_Array();
    PackedArray2D_Array(quint32 sz);
    PackedArray2D_Array(quint32 sz, const T &value);
    PackedArray2D_Array(const QVector<T> &values);
    
    PackedArray2D_Array(const PackedArray2D_Array<T> &other);
    
    ~PackedArray2D_Array();
    
    PackedArray2D_Array<T>& operator=(const PackedArray2D_Array<T> &other);
    
    bool operator==(const PackedArray2D_Array<T> &other) const;
    bool operator!=(const PackedArray2D_Array<T> &other) const;

    const T& operator[](quint32 i) const;
    T& operator[](quint32 i);
    
    const T& at(quint32 i) const;

    int count() const;
    int size() const;
    
    int nValues() const;
    
    bool isEmpty() const;
    
    const T* data() const;
    T* data();
    
    const T* constData() const;
    
    void update(const PackedArray2D_Array<T> &other);
    
    void assertValidIndex(quint32 i) const;

protected:
    PackedArray2D_Array(detail::PackedArray2D_ArrayData<T> *data);
    
private:
    /** Implicitly shared pointer to the array */
    detail::SharedArray2DPtr< detail::PackedArray2D_ArrayData<T> > d;
};

}

/** This class provides an array of arrays of type T, where
    each array is packed in memory so that each object of type
    T is contiguous. This allows rapid indexing over each
    object in each array, or over all objects in all arrays.
    
    @author Christopher Woods
*/
template<class T>
class PackedArray2D
{

friend class detail::PackedArray2DMemory<T>;

public:
    typedef typename detail::PackedArray2D_Array<T> Array;

    PackedArray2D();
    
    PackedArray2D(const Array &array);
    PackedArray2D(const QVector<Array> &arrays);
    
    PackedArray2D(const QVector<T> &values);
    PackedArray2D(const QVector< QVector<T> > &values);
    
    PackedArray2D(const PackedArray2D<T> &other);
    
    ~PackedArray2D();
    
    PackedArray2D<T>& operator=(const PackedArray2D<T> &other);
    
    bool operator==(const PackedArray2D<T> &other) const;
    bool operator!=(const PackedArray2D<T> &other) const;
    
    const Array& operator[](quint32 i) const;
    
    const T& operator()(quint32 i, quint32 j) const;
    T& operator()(quint32 i, quint32 j);
    
    const Array& at(quint32 i) const;
    const T& at(quint32 i, quint32 j) const;
    
    int count() const;
    int size() const;
    
    int nArrays() const;
    int nValues() const;

    bool isEmpty() const;

    const Array* data() const;
    const Array* constData() const;
    
    const T* data(quint32 i) const;
    T* data(quint32 i);
    
    const T* constData(quint32 i) const;
    
    const T* valueData() const;
    T* valueData();
    
    const T* constValueData() const;

    void update(quint32 i, const Array &array);

    void assertValidIndex(quint32 i) const;

private:
    /** Implicitly shared pointer to the array data */
    detail::SharedArray2DPtr< detail::PackedArray2DData<T> > d;
};

}

namespace SireBase
{

namespace detail
{

////////
//////// Implementation of PackedArray2DMemory
////////

/** Pointer to the global null PackedArray */
template<class T>
SharedArray2DPtr< PackedArray2DData<T> > PackedArray2DMemory<T>::shared_null;

/** Create space for the arrays. The layout in memory is;

    --------------------------------------------------
    | A |   B   |    C    |            D             |
    --------------------------------------------------

      A = PackedArray2DData<T> object          1 * sizeof(PackedArray2DData<T>)
      B = PackedArray2D_Array<T> objects      narrays * sizeof(PackedArray2D_Array<T>)
      C = PackedArray2D_ArrayData<T> objects   narrays * sizeof( " _ArrayData<T>)
      D = Array of T objects                 nvalues * sizeof(T)
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
char* PackedArray2DMemory<T>::create(quint32 narrays, quint32 nvalues)
{
    //create space for the data
    char *storage = PackedArray2DMemoryBase::create(narrays, nvalues, 
                                            sizeof(PackedArray2DData<T>),
                                            sizeof(PackedArray2D_Array<T>),
                                            sizeof(PackedArray2D_ArrayData<T>),
                                            sizeof(T));
    
    try
    {
        quint32 sz = PackedArray2DMemoryBase::getSize(narrays, nvalues,   
                                            sizeof(PackedArray2DData<T>),
                                            sizeof(PackedArray2D_Array<T>),
                                            sizeof(PackedArray2D_ArrayData<T>),
                                            sizeof(T));
    
        //now we need to construct each object in turn at its correct
        //location in this space
        PackedArray2DData<T> *arraydata 
                    = new (storage) PackedArray2DData<T>(narrays, nvalues);
                    
        //add a reference count for this object
        arraydata->ref.ref();

        //advance the index into the storage array so that we are now
        //pointing just after the PackedArray2DData
        quint32 idx = sizeof(PackedArray2DData<T>);
                    
        //this is the location of the first PackedArray2D_Array<T> :-)
        PackedArray2DMemoryBase::setArray0(arraydata, idx);
        quint32 array0_idx = idx;
        
        //the first PackedArray2D_ArrayData lies at this 
        //index + narrays*sizeof(PackedArray2D_Array<T>)
        quint32 dataidx = idx + narrays * sizeof(PackedArray2D_Array<T>);
        
        PackedArray2DMemoryBase::setArrayData0(arraydata, dataidx);

        if (narrays > 0)
        {
            //loop over each array and create it in its place
            for (quint32 i=0; i<narrays; ++i)
            {   
                //assert that there is sufficient space in the array
                BOOST_ASSERT(idx + sizeof(PackedArray2D_Array<T>) < sz);
                BOOST_ASSERT(dataidx + sizeof(PackedArray2D_ArrayData<T>) < sz);
            
                //create the ArrayData, letting it know where it
                //is relative to the beginning of the storage array
                PackedArray2D_ArrayData<T> *data 
                        = new (storage + dataidx) PackedArray2D_ArrayData<T>(dataidx);
                
                //now create the PackedArray2D_Array<T> that uses this data
                new (storage + idx) PackedArray2D_Array<T>(data);
        
                //advance the index into the storage array to point
                //just after the just-created PackedArray2D_ArrayData
                idx += sizeof(PackedArray2D_Array<T>);
                dataidx += sizeof(PackedArray2D_ArrayData<T>);
            }
            
            //make idx point just after all of the PackedArray2D_ArrayDatas...
            BOOST_ASSERT( idx == array0_idx 
                                + narrays*sizeof(PackedArray2D_Array<T>) );
                                
            idx = dataidx;
        }
        else
        {
            //we need to create space for the null PackedArray2D_Array<T>
            quint32 data_idx = idx + sizeof(PackedArray2D_Array<T>);
            PackedArray2D_ArrayData<T> *array_arraydata 
                        = new (storage + data_idx) PackedArray2D_ArrayData<T>(data_idx);
            
            new (storage + idx) PackedArray2D_Array<T>(array_arraydata);
            
            PackedArray2DMemoryBase::setNValues(array_arraydata, 0);
            PackedArray2DMemoryBase::setValue0(array_arraydata, 0);
            
            idx += sizeof(PackedArray2D_Array<T>) 
                    + sizeof(PackedArray2D_ArrayData<T>);
        }
        
        //we are now at the location of the first item in the array
        PackedArray2DMemoryBase::setValue0(arraydata, idx);
        
        //loop over each object and create it in place
        for (quint32 i=0; i<nvalues; ++i)
        {
            BOOST_ASSERT(idx + sizeof(T) <= sz);

            new (storage + idx) T();
                
            idx += sizeof(T);
        }
        
        //we should now be at the end of the storage
        BOOST_ASSERT( idx == sz );
                
        return storage;
    }
    catch(...)
    {
        delete[] storage;
        throw;
    }
}

/** Increase the reference count */
template<class T>
SIRE_INLINE_TEMPLATE
void PackedArray2DMemory<T>::incref(char *this_ptr, quint32 this_idx)
{
    ( (PackedArray2DData<T>*)(PackedArray2DMemoryBase::getRoot(this_ptr,this_idx)) )
                    ->incref();
}

/** Decrease the reference count */
template<class T>
SIRE_INLINE_TEMPLATE
void PackedArray2DMemory<T>::decref(char *this_ptr, quint32 this_idx)
{
    ( (PackedArray2DData<T>*)(PackedArray2DMemoryBase::getRoot(this_ptr,this_idx)) )
                    ->decref();
}

/** Destroy the object 'array' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2DMemory<T>::destroy(PackedArray2DData<T> *array)
{
    //we need to delete it in the opposite order to creation - first
    //lets delete all of the objects
    quint32 nvalues = array->nValues();

    char *storage = (char*)array;
        
    if (nvalues > 0)
    {
        T *values = array->valueData();
        
        for (qint32 i=nvalues-1; i>=0; --i)
        {
            values[i].~T();
        }
    }
        
    //now delete all of the arrays
    quint32 narrays = array->nArrays();
        
    if (narrays > 0)
    {
        PackedArray2D_Array<T> *arrays = array->arrayData();
            
        for (qint32 i=narrays-1; i>=0; --i)
        {
            PackedArray2D_ArrayData<T> *arraydata 
                 = const_cast<PackedArray2D_ArrayData<T>*>( arrays[i].d.constData() );
            
            //delete the PackedArray2D_ArrayData
            arraydata->~PackedArray2D_ArrayData<T>();
            
            //remove the shared pointer
            arrays[i].d.weakRelease();
            
            //delete the PackedArray2D_Array<T>
            arrays[i].~PackedArray2D_Array<T>();
        }
    }
        
    //delete the PackedArray2DData object itself
    array->~PackedArray2DData<T>();
        
    //finally release the memory
    delete[] storage;
}

/** Detach the data pointed to by 'this_ptr' from shared storage */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
char* PackedArray2DMemory<T>::detach(char *this_ptr, quint32 this_idx)
{
    //get a pointer to the start of the storage for this container
    char *storage = getRoot(this_ptr, this_idx);
    
    //The PackedArray2DData object is at the beginning of this storage array
    PackedArray2DData<T> *arraydata = (PackedArray2DData<T>*) storage;
    
    if (arraydata->ref != 1)
    {
        //there is more than one reference to this data - it will have to 
        //be cloned - get the size of memory to be cloned
        int sz = getSize(arraydata->nArrays(), arraydata->nValues(),
                         sizeof(PackedArray2DData<T>), 
                         sizeof(PackedArray2D_Array<T>),
                         sizeof(PackedArray2D_ArrayData<T>),
                         sizeof(T));
                         
        //create space for the copy
        char *new_storage = new char[sz];
        
        //copy the data
        void *output = qMemCopy( new_storage, storage, sz*sizeof(char) );
        BOOST_ASSERT(output == new_storage);

        //the first part of the data is the PackedArray2DData object
        PackedArray2DData<T> *new_arraydata = (PackedArray2DData<T>*) new_storage;
        
        //set the reference count of this copy to 1
        new_arraydata->ref = 1;

        //now loose a reference to the original
        PackedArray2DMemory<T>::decref(this_ptr, this_idx);
        
        //the final step is to update all of the PackedArray2D_Array<T>
        //pointers that exist in this array (otherwise they will all
        //point to the original!)
        
        if (new_arraydata->nArrays() > 0)
        {
            PackedArray2D_Array<T> *arrays = new_arraydata->arrayData();
            PackedArray2D_ArrayData<T> *arrays_data = new_arraydata->arrayDataData();
        
            for (quint32 i=0; i < new_arraydata->nArrays(); ++i)
            {
                arrays[i].d.weakRelease();
                arrays[i].d.weakAssign( &(arrays_data[i]) );
            }
        }
        else
        {
            //update the null PackedArray2D_Array<T>
            PackedArray2D_Array<T> *arrays = new_arraydata->arrayData();
            PackedArray2D_ArrayData<T> *arrays_data = new_arraydata->arrayDataData();

            arrays[0].d.weakRelease();
            arrays[0].d.weakAssign( &(arrays_data[0]) );
        }

        //return a pointer to the clone
        return new_storage + this_idx;
    }
    else
    {
        //only one reference, so no need to clone
        return this_ptr;
    }
}

///////
/////// Implementation of PackedArray2DData<T>
///////

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2DData<T>::PackedArray2DData() : PackedArray2DDataBase()
{}

/** Construct to hold the specified number of arrays and values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2DData<T>::PackedArray2DData(quint32 narrays, quint32 nvalues)
                   : PackedArray2DDataBase(narrays, nvalues)
{}

/** Copy constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2DData<T>::PackedArray2DData(const PackedArray2DData &other)
                   : PackedArray2DDataBase(other)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2DData<T>::~PackedArray2DData()
{}

/** Decrease the reference count of this object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2DData<T>::decref()
{
    if (not this->ref.deref())
        PackedArray2DMemory<T>::destroy(this);
}

/** Return a pointer to a copy of this array that has been 
    detached from shared storage */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2DData<T>* PackedArray2DData<T>::detach()
{
    return (PackedArray2DData<T>*)( PackedArray2DMemory<T>::detach( (char*)(this), 0 ) );
}

/** Return a pointer to the null array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const PackedArray2D_ArrayData<T>* PackedArray2DData<T>::nullArray() const
{
    BOOST_ASSERT( this->nArrays() == 0 );
    BOOST_ASSERT( this->getArrayData0() != 0 );
    
    return (const PackedArray2D_ArrayData<T>*)
                    ( this->memory() + this->getArrayData0() );
}

/** Return a raw pointer to the array of PackedArray2D_ArrayData objects */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const PackedArray2D_ArrayData<T>* PackedArray2DData<T>::arrayDataData() const
{
    if (this->getArrayData0() == 0)
        return 0;
    else
        return (const PackedArray2D_ArrayData<T>*)
                    ( this->memory() + this->getArrayData0() );
}

/** Return a raw pointer to the array of PackedArray2D_Array<T> objects */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const PackedArray2D_Array<T>* PackedArray2DData<T>::arrayData() const
{
    if (this->getArray0() == 0)
        return 0;
    else
        return (const PackedArray2D_Array<T>*)( this->memory() + this->getArray0() );
}

/** Return a raw pointer to the array of objects in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PackedArray2DData<T>::valueData() const
{
    if (this->getValue0() == 0)
        return 0;
    else
        return (const T*)( this->memory() + this->getValue0() );
}

/** Return a raw pointer to the array of PackedArray2D_ArrayData objects */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_ArrayData<T>* PackedArray2DData<T>::arrayDataData()
{
    if (this->getArrayData0() == 0)
        return 0;
    else
        return (PackedArray2D_ArrayData<T>*)( this->memory() + this->getArrayData0() );
}

/** Return a raw pointer to the array of PackedArray2D_Array<T> objects */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_Array<T>* PackedArray2DData<T>::arrayData()
{
    if (this->getArray0() == 0)
        return 0;
    else
        return (PackedArray2D_Array<T>*)( this->memory() + this->getArray0() );
}

/** Return a raw pointer to the array of objects in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T* PackedArray2DData<T>::valueData()
{
    if (this->getValue0() == 0)
        return 0;
    else
        return (T*)( this->memory() + this->getValue0() );
}

/** Set the number of values in the ith array - can only be called
    before the close() (so don't call it after closing the arrays) */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2DData<T>::setNValuesInArray(quint32 i, quint32 nvalues)
{
    BOOST_ASSERT( i < this->nArrays() );
    PackedArray2DMemoryBase::setNValues(&(this->arrayDataData()[i]), nvalues);
}

/** Close out this array - this fixes the number of values
    in each array in this array of arrays */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2DData<T>::close()
{
    PackedArray2D_ArrayData<T> *arraydata = this->arrayDataData();
    
    //ensure that everything adds up
    quint32 n_assigned_vals = 0;
    
    quint32 value0 = this->getValue0();
    
    //close out each array
    for (quint32 i=0; i<this->nArrays(); ++i)
    {
        PackedArray2D_ArrayData<T> &array = arraydata[i];
        
        //tell the array where its first value is located
        if (array.nValues() > 0)
            PackedArray2DMemoryBase::setValue0(&array, 
                                        value0 + n_assigned_vals * sizeof(T));
        else
            PackedArray2DMemoryBase::setValue0(&array, 0);

        n_assigned_vals += array.nValues();
    }

    BOOST_ASSERT( n_assigned_vals == this->nValues() );
}

////////
//////// Implementation of PackedArray2D_ArrayData
////////

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_ArrayData<T>::PackedArray2D_ArrayData()
                         : PackedArray2D_ArrayDataBase()
{}

/** Construct the data, telling it that it is at 
    index 'this_idx' in the memory storage array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_ArrayData<T>::PackedArray2D_ArrayData(quint32 this_idx)
                         : PackedArray2D_ArrayDataBase(this_idx)
{}

/** Copy constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_ArrayData<T>::PackedArray2D_ArrayData(const PackedArray2D_ArrayData<T> &other)
                         : PackedArray2D_ArrayDataBase(other)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_ArrayData<T>::~PackedArray2D_ArrayData()
{}

/** Return a pointer to a PackedArray2D_ArrayData that has been detached from 
    shared storage (could be a pointer to this PackedArray2D_ArrayData) */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_ArrayData<T>* PackedArray2D_ArrayData<T>::detach()
{
    return (PackedArray2D_ArrayData<T>*)
           ( PackedArray2DMemory<T>::detach( (char*)this, this->getThisArray() ) );
}

/** Return a pointer to a PackedArray2DData object that contains 
    just the data of this PackedArray2D_Array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2DData<T>* PackedArray2D_ArrayData<T>::extract() const
{
    const PackedArray2DData<T> *array = (const PackedArray2DData<T>*)(this->memory());
    
    if (array->nArrays() <= 1)
        //this is already extracted!
        return const_cast<PackedArray2DData<T>*>(array);
    
    char *new_storage = PackedArray2DMemory<T>::create(1, this->nValues());
    
    try
    {
                
    //ok, we need to extract!
    PackedArray2DData<T> *new_array = (PackedArray2DData<T>*)(new_storage);

    new_array->setNValuesInArray(1, this->nValues());
    new_array->close();
    
    //copy the objects
    void *output = qMemCopy( new_array->valueData(),
                             this->valueData(),
                             this->nValues() * sizeof(T) );
                       
    BOOST_ASSERT( output == new_array->valueData() );

    return new_array;

    }
    catch(...)
    {
        delete[] new_storage;
    }
    
    return 0;
}

/** Increment the reference count for this object */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D_ArrayData<T>::incref()
{
    PackedArray2DMemory<T>::incref( (char*)this, this->getThisArray() );
}

/** Decrement the reference count for this object - this may delete it! */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D_ArrayData<T>::decref()
{
    PackedArray2DMemory<T>::decref( (char*)this, this->getThisArray() );
}

/** Return a pointer to the first object in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PackedArray2D_ArrayData<T>::valueData() const
{
    if (this->getValue0() == 0)
        return 0;
    else
        return (const T*)( this->memory() + this->getValue0() );
}

/** Return a pointer to the first object in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T* PackedArray2D_ArrayData<T>::valueData()
{
    if (this->getValue0() == 0)
        return 0;
    else
        return (T*)( this->memory() + this->getValue0() );
}

///////////
/////////// Implementation of PackedArray2D_Array
///////////

template<class T>
static const SharedArray2DPtr< PackedArray2DData<T> >& getSharedNull()
{
    if (detail::PackedArray2DMemory<T>::shared_null.constData() == 0)
    {
        detail::PackedArray2DMemory<T>::shared_null 
                = (PackedArray2DData<T>*)( PackedArray2DMemory<T>::create(0,0) );
                
        detail::PackedArray2DMemory<T>::shared_null->close();
    }
    
    return detail::PackedArray2DMemory<T>::shared_null;
}

template<class T>
static SharedArray2DPtr< PackedArray2D_ArrayData<T> > createNullArray()
{
    const SharedArray2DPtr< PackedArray2DData<T> > &array 
                            = SireBase::detail::getSharedNull<T>();
    
    return SharedArray2DPtr< PackedArray2D_ArrayData<T> >( array->nullArray() );
}

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_Array<T>::PackedArray2D_Array()
                     : d( SireBase::detail::createNullArray<T>() )
{}

template<class T>
static SharedArray2DPtr< PackedArray2D_ArrayData<T> > createArray(quint32 sz)
{
    if (sz == 0)
        return SireBase::detail::createNullArray<T>();
        
    //construct space for 1 array of sz objects
    char *storage = PackedArray2DMemory<T>::create(1, sz);
        
    PackedArray2DData<T> *array = (PackedArray2DData<T>*)storage;
        
    array->setNValuesInArray(0, sz);
    array->close();
    
    return SharedArray2DPtr< PackedArray2D_ArrayData<T> >( &(array->arrayDataData()[0]) );
}

/** Construct an array of size 'sz'. All of the objects will be initialised  
    with their default values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_Array<T>::PackedArray2D_Array(quint32 sz)
                     : d( SireBase::detail::createArray<T>(sz) )
{}

/** Construct an array of size 'sz' where all of the objects have value 'value'. */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_Array<T>::PackedArray2D_Array(quint32 sz, const T &value)
                     : d( SireBase::detail::createArray<T>(sz) )
{
    if (sz == 0)
        return;

    T *data = d->valueData();

    //initialise all of the values
    for (quint32 i=0; i<sz; ++i)
    {
        data[i] = value;
    }
}

/** Construct from the passed array of values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_Array<T>::PackedArray2D_Array(const QVector<T> &values)
{
    PackedArray2D<T> array(values);
    this->operator=( array.constData()[0] );
}

/** Private constructor used internally */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_Array<T>::PackedArray2D_Array(detail::PackedArray2D_ArrayData<T> *data)
{
    d.weakAssign(data);
}

/** Copy constructor - fast as this class is implicitly shared */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_Array<T>::PackedArray2D_Array(const PackedArray2D_Array<T> &other)
                     : d(other.d)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_Array<T>::~PackedArray2D_Array()
{}

/** Copy assignment operator - fast as this class is implicitly shared.
    However, as this copies a view of the whole packed array, you cannot
    use this to assign the contents of this array within its current 
    PackedArray2D - to do this, you must use the 'update()' function. */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D_Array<T>& PackedArray2D_Array<T>::operator=(
                                        const PackedArray2D_Array<T> &other)
{
    this->d = other.d;
    return *this;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PackedArray2D_Array<T>::operator==(const PackedArray2D_Array<T> &other) const
{
    if (this->d.constData() == other.d.constData())
        return true;

    if (this->nValues() != other.nValues())
        return false;
        
    const T *this_data = this->constData();
    const T *other_data = other.constData();
    
    for (quint32 i=0; i<this->nValues(); ++i)
    {
        if (this_data[i] != other_data[i])
            return false;
    }
    
    return true;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PackedArray2D_Array<T>::operator!=(const PackedArray2D_Array<T> &other) const
{
    return not this->operator==(other);
}

/** Assert that the index 'i' is valid for this array 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D_Array<T>::assertValidIndex(quint32 i) const
{
    if (i >= this->nValues())
        detail::throwPackedArray2D_Array_invalidIndex(i, this->nValues());
}

/** Return a raw pointer to the array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PackedArray2D_Array<T>::data() const
{
    return this->d->valueData();
}

/** Return a modifiable raw pointer to the array
     *DO NOT* delete this pointer! */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T* PackedArray2D_Array<T>::data()
{
    return this->d->valueData();
}

/** Return a raw pointer to the array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PackedArray2D_Array<T>::constData() const
{
    return this->d->valueData();
}

/** Return a reference to the ith object in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& PackedArray2D_Array<T>::operator[](quint32 i) const
{
    this->assertValidIndex(i);
    return this->data()[i];
}

/** Return a reference to the ith object in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T& PackedArray2D_Array<T>::operator[](quint32 i)
{
    this->assertValidIndex(i);
    return this->data()[i];
}

/** Return a reference to the ith object in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& PackedArray2D_Array<T>::at(quint32 i) const
{
    this->assertValidIndex(i);
    return this->data()[i];
}

/** Return the number of objects in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int PackedArray2D_Array<T>::count() const
{
    return d->nValues();
}

/** Return the number of objects in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int PackedArray2D_Array<T>::size() const
{
    return this->d->nValues();
}

/** Return the number of objects in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int PackedArray2D_Array<T>::nValues() const
{
    return this->d->nValues();
}

/** Return whether or not this array is empty */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PackedArray2D_Array<T>::isEmpty() const
{
    return this->count() == 0;
}

/** Update this array so that it has the same contents as 'other'

    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D_Array<T>::update(const PackedArray2D_Array<T> &other)
{
    if (this->size() != other.size())
        detail::throwPackedArray2D_Array_incompatibleError(this->size(),other.size());
        
    if (this->size() == 0)
        return;
        
    void *output = qMemCopy(this->data(), other.data(), 
                            this->size() * sizeof(T));
                            
    BOOST_ASSERT( output == this->data() );
}

} // end of namespace detail

/////////
///////// Implementation of PackedArray2D
/////////

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<T>::PackedArray2D() : d( SireBase::detail::getSharedNull<T>() )
{}

/** Construct from a passed Array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<T>::PackedArray2D(const PackedArray2D<T>::Array &array)
               : d( array.d->extract() )
{}

/** Copy assignment operator - fast as this class is implicitly shared */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<T>& PackedArray2D<T>::operator=(const PackedArray2D<T> &other)
{
    d = other.d;
    return *this;
}

namespace detail
{

template<class T>
static SharedArray2DPtr< PackedArray2DData<T> > 
createArray(quint32 narrays, quint32 nvals)
{
    if (narrays == 0)
        return SireBase::detail::getSharedNull<T>();
        
    //construct space for narrays arrays of nvals objects
    char *storage = PackedArray2DMemory<T>::create(narrays, nvals);
        
    PackedArray2DData<T> *array = (PackedArray2DData<T>*)storage;
    
    return SharedArray2DPtr< PackedArray2DData<T> >(array);
}

}

/** Construct from an array of arrays */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<T>::PackedArray2D(const QVector<PackedArray2D<T>::Array> &arrays)
{
    if (arrays.isEmpty())
        return;

    else if (arrays.count() == 1)
    {
        this->operator=( PackedArray2D<T>(arrays.at(0)) );
        return;
    }

    //count the number of values and groups...
    quint32 nvals = 0;
    
    const Array *arrays_data = arrays.constData();
    quint32 narrays = arrays.count();
    
    for (quint32 i=0; i<narrays; ++i)
    {
        nvals += arrays_data[i].count();
    }
    
    d = SireBase::detail::createArray<T>(narrays, nvals);
    
    detail::PackedArray2DData<T> *dptr = d.data();
    T *values = dptr->valueData();
    
    //dimension each packed array
    for (quint32 i=0; i<narrays; ++i)
    {
        //dimension the array
        dptr->setNValuesInArray(i, arrays_data[i].count());
    }
    
    dptr->close();
    
    //now copy all of the data
    for (quint32 i=0; i<narrays; ++i)
    {
        const Array &array = arrays_data[i];
        const T *array_values = array.constData();
        
        void *output = qMemCopy(values, array_values, array.count()*sizeof(T));
        
        BOOST_ASSERT( output == values );
        
        values += array.count();
    }
}

/** Construct from an array of values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<T>::PackedArray2D(const QVector<T> &values)
{
    this->operator=( PackedArray2D<T>(Array(values)) ); 
}

/** Construct from an array of array of values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<T>::PackedArray2D(const QVector< QVector<T> > &values)
{
    if (values.isEmpty())
        return;
        
    else if (values.count() == 1)
    {
        this->operator=( PackedArray2D<T>(values.at(0)) );
        return;
    }

    //count the number of values and groups...
    quint32 nvals = 0;
    
    const QVector<T> *values_data = values.constData();
    quint32 narrays = values.count();
    
    for (quint32 i=0; i<narrays; ++i)
    {
        nvals += values_data[i].count();
    }
    
    d = SireBase::detail::createArray<T>(narrays, nvals);
    
    detail::PackedArray2DData<T> *dptr = d.data();
    
    //dimension each packed array
    for (quint32 i=0; i<narrays; ++i)
    {
        //dimension the array
        dptr->setNValuesInArray(i, values_data[i].count());
    }
    
    dptr->close();
    
    //now copy all of the data
    T *values_array = dptr->valueData();

    for (quint32 i=0; i<narrays; ++i)
    {
        const QVector<T> &array = values_data[i];
        const T *array_values = array.constData();
        
        void *output = qMemCopy(values_array, array_values, array.count()*sizeof(T));
        
        BOOST_ASSERT( output == values_array );
        
        values_array += array.count();
    }
}

/** Copy constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<T>::PackedArray2D(const PackedArray2D<T> &other)
               : d(other.d)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
PackedArray2D<T>::~PackedArray2D()
{}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PackedArray2D<T>::operator==(const PackedArray2D<T> &other) const
{
    if (d.constData() == other.d.constData())
        return true;
        
    if (this->nArrays() != other.nArrays())
        return false;
        
    const Array *this_arrays = this->constData();
    const Array *other_arrays = other.constData();
    quint32 narrays = this->count();
    
    for (quint32 i=0; i<narrays; ++i)
    {
        if (this_arrays[i] != other_arrays[i])
            return false;
    }
    
    return true;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PackedArray2D<T>::operator!=(const PackedArray2D<T> &other) const
{
    return not this->operator==(other);
}

/** Assert that the index 'i' is valid for this array

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D<T>::assertValidIndex(quint32 i) const
{
    if (i >= this->count())
        SireBase::detail::throwPackedArray2D_invalidIndex(i, this->count());
}

/** Return a reference to the ith array in this container

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename PackedArray2D<T>::Array& PackedArray2D<T>::operator[](quint32 i) const
{
    this->assertValidIndex(i);
    return d->arrayData()[i];
}

/** Return a reference to the ith array in this container

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename PackedArray2D<T>::Array& PackedArray2D<T>::at(quint32 i) const
{
    return this->operator[](i);
}

/** Return a reference to the jth value of the ith array in 
    this container
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& PackedArray2D<T>::operator()(quint32 i, quint32 j) const
{
    return this->at(i).at(j);
}

/** Return a modifiable reference to the jth value of the ith array 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T& PackedArray2D<T>::operator()(quint32 i, quint32 j)
{
    this->assertValidIndex(i);
    Array &array = d->arrayData()[i];
    
    return array[j];
}

/** Return a reference to the jth value of the ith array in 
    this container
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& PackedArray2D<T>::at(quint32 i, quint32 j) const
{
    return this->operator()(i,j);
}

/** Return the number of arrays in this PackedArray2D */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int PackedArray2D<T>::count() const
{
    return d->nArrays();
}

/** Return the number of arrays in this PackedArray2D */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int PackedArray2D<T>::size() const
{
    return this->count();
}

/** Return the number of arrays in this PackedArray2D */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int PackedArray2D<T>::nArrays() const
{
    return this->count();
}

/** Return the total number of objects in this array */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int PackedArray2D<T>::nValues() const
{
    return d->nValues();
}

/** Return whether or not this array is empty (contains no arrays) */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool PackedArray2D<T>::isEmpty() const
{
    return this->count() == 0;
}

/** Return a raw pointer to the array of arrays. */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename PackedArray2D<T>::Array* PackedArray2D<T>::data() const
{
    return d->arrayData();
}

/** Return a raw pointer to the array of arrays. */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename PackedArray2D<T>::Array* PackedArray2D<T>::constData() const
{
    return d->arrayData();
}

/** Return a raw pointer to the array of values in the ith array 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PackedArray2D<T>::data(quint32 i) const
{
    return this->at(i).data();
}

/** Return a modifiable pointer to the array of values in the
    ith array. *DO NOT* delete this pointer! 
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T* PackedArray2D<T>::data(quint32 i)
{
    this->assertValidIndex(i);
    return d->arrayData()[i].data();
}

/** Return a raw pointer to the array of values in the ith array 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PackedArray2D<T>::constData(quint32 i) const
{
    return this->data(i);
}

/** Return a raw pointer to the array of all values in all arrays */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PackedArray2D<T>::valueData() const
{
    return d->valueData();
}

/** Return a modifiable pointer to the array of all values in 
    all arrays. *DO NOT* delete this pointer. */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T* PackedArray2D<T>::valueData()
{
    return d->valueData();
}

/** Return a raw pointer to the array of all values in all arrays */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* PackedArray2D<T>::constValueData() const
{
    return d->valueData();
}

/** Update the ith array so that it has the same contents as 'array'

    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void PackedArray2D<T>::update(quint32 i, const Array &array)
{
    this->assertValidIndex(i);
    d->arrayData()[i].update(array);
}

} // end of namespace SireBase

SIRE_END_HEADER

#endif
