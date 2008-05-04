/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#ifndef SIREMOL_ATOMPROPERTY_HPP
#define SIREMOL_ATOMPROPERTY_HPP

#include <QVector>
#include <QVariant>

#include "molviewproperty.h"
#include "moleculeinfodata.h"
#include "atomselection.h"

#include "SireBase/packedarray2d.hpp"

#include "SireError/errors.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class AtomProp;

template<class T>
class AtomProperty;
}

QDataStream& operator<<(QDataStream&, const SireMol::AtomProp&);
QDataStream& operator>>(QDataStream&, SireMol::AtomProp&);

template<class T>
QDataStream& operator<<(QDataStream&, const SireMol::AtomProperty<T>&);
template<class T>
QDataStream& operator>>(QDataStream&, SireMol::AtomProperty<T>&);

namespace SireMol
{

using SireBase::Property;
using SireBase::PackedArray2D;

/** Small class used to give a common base to all
    AtomProperty classes */
class SIREMOL_EXPORT AtomProp : public MolViewProperty
{
public:
    AtomProp();
    AtomProp(const AtomProp &other);
    
    virtual ~AtomProp();
    
    virtual bool canConvert(const QVariant &value) const=0;
    
    virtual void assignFrom(const QVector<QVariant> &values)=0;
    virtual void assignFrom(const QVector< QVector<QVariant> > &values)=0;
    
    virtual QVector< QVector<QVariant> > toVariant() const=0;
    
    virtual void assertCanConvert(const QVariant &value) const=0;
};

/** This is a property that can hold one value for each
    atom in the molecule. The properties are held in 
    a packed array, with each array corresponding to a CutGroup,
    and the order of the properties stored according 
    to CutGroup
    
    mol.setProperty( "charge", AtomCharges( [....] ) )
    mol.setProperty( "lj", AtomLJs( [ [...], [...] ] ) )

    atom.setProperty( "charge", 0.0 * mod_electron )
    
    @author Christopher Woods
*/
template<class T>
class AtomProperty 
    : public SireBase::ConcreteProperty<AtomProperty<T>, AtomProp>
{

friend QDataStream& ::operator<<<>(QDataStream&, const AtomProperty<T>&);
friend QDataStream& ::operator>><>(QDataStream&, AtomProperty<T>&);

public:
    AtomProperty();

    AtomProperty(const MoleculeInfoData &molinfo);
    
    AtomProperty(const T &value);
    AtomProperty(const PackedArray2D<T> &values);

    AtomProperty(const QVector< QVector<QVariant> > &values);
    AtomProperty(const QVector<QVariant> &values);
    
    AtomProperty(const AtomProperty<T> &other);
    
    ~AtomProperty();
    
    AtomProperty<T>& operator=(const AtomProperty<T> &other);
    
    static const char* typeName()
    {
        return QMetaType::typeName( qMetaTypeId< AtomProperty<T> >() );
    }
    
    AtomProperty<T>* clone() const
    {
        return new AtomProperty<T>(*this);
    }
 
    bool operator==(const AtomProperty<T> &other) const;
    bool operator!=(const AtomProperty<T> &other) const;

    QVector< QVector<QVariant> > toVariant() const;
    
    void assignFrom(const QVector<QVariant> &values);
    void assignFrom(const QVector< QVector<QVariant> > &values);

    const typename PackedArray2D<T>::Array& operator[](CGIdx cgidx) const;
    const typename PackedArray2D<T>::Array& at(CGIdx cgidx) const;
    const typename PackedArray2D<T>::Array& get(CGIdx cgidx) const;

    const T& operator[](const CGAtomIdx &cgatomidx) const;
    const T& at(const CGAtomIdx &cgatomidx) const;
    const T& get(const CGAtomIdx &cgatomidx) const;

    AtomProperty<T>& set(const CGAtomIdx &cgatomidx, const T &value);

    AtomProperty<T>& set(CGIdx cgidx, const QVector<T> &values);

    const PackedArray2D<T>& array() const;

    const typename PackedArray2D<T>::Array* data() const;
    const typename PackedArray2D<T>::Array* constData() const;

    const T* data(CGIdx cgidx) const;
    const T* constData(CGIdx cgidx) const;

    int size() const;
    int count() const;
    
    int nCutGroups() const;
    
    int nAtoms() const;
    int nAtoms(CGIdx cgidx) const;

    bool isCompatibleWith(const MoleculeInfoData &molinfo) const;

    AtomProperty<T> matchToSelection(const AtomSelection &selection) const;

    bool canConvert(const QVariant &value) const;

    void assertCanConvert(const QVariant &value) const;

private:
    /** The actual atomic property values */
    SireBase::PackedArray2D<T> props;
};

/** Return whether or not the variant 'value' can be converted to be
    held as an AtomProperty<T> */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool AtomProperty<T>::canConvert(const QVariant &value) const
{
    return value.canConvert<T>();
}

/** Assert that the passed variant value can be converted to be
    held within this property
    
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomProperty<T>::assertCanConvert(const QVariant &value) const
{
    if (not value.canConvert<T>())
    {
        throw SireError::invalid_cast( QObject::tr(
            "It is not possible to convert the value of type %1 to "
            "type %2, as is required for storing in the AtomProperty %3.")
                .arg(value.typeName()).arg(T::typeName())
                .arg(this->what()), CODELOC );
    }
}

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::AtomProperty()
                : SireBase::ConcreteProperty<AtomProperty<T>,AtomProp>()
{}

/** Create an AtomProperty that holds one value for each 
    atom described in 'molinfo'. Each atom starts with
    a default-constructed value of the property */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::AtomProperty(const MoleculeInfoData &molinfo)
                : SireBase::ConcreteProperty<AtomProperty<T>,AtomProp>()
{   
    int ncg = molinfo.nCutGroups();

    if (ncg > 0)
    {
        //create space for each CutGroup
        QVector< QVector<T> > tmp_props = QVector< QVector<T> >(ncg);
        
        QVector<T> *tmp_props_array = tmp_props.data();
        
        for (CGIdx i(0); i<ncg; ++i)
        {
            //now create space for all of the atoms
            tmp_props_array[i] = QVector<T>(molinfo.nAtoms(i));
        }

        //now copy this into the PackedArray
        props = PackedArray2D<T>(tmp_props);
    }
}

/** Construct an Atom property that holds a single value (only
    suitable for a molecule that has just one atom in just one CutGroup) */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::AtomProperty(const T &value)
                : SireBase::ConcreteProperty<AtomProperty<T>,AtomProp>()
{
    QVector<T> tmp_props(1, value);
    props = PackedArray2D<T>(tmp_props);
}

/** Construct the Atom property from the PackedArray2D of values */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::AtomProperty(const PackedArray2D<T> &values)
                : SireBase::ConcreteProperty<AtomProperty<T>,AtomProp>(),
                  props(values)
{}

/** Construct the Atom property from the array of QVariants.
    
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::AtomProperty(const QVector<QVariant> &values)
                : SireBase::ConcreteProperty<AtomProperty<T>,AtomProp>()
{
    if (values.isEmpty())
        return;
    
    int nvals = values.count();
    const QVariant *values_array = values.constData();
    
    QVector<T> converted_vals(nvals);
    T *converted_vals_array = converted_vals.data();
    
    for (int i=0; i<nvals; ++i)
    {
        const QVariant &value = values_array[i];
        AtomProperty<T>::assertCanConvert(value);

        converted_vals_array[i] = value.value<T>();
    }
    
    props = PackedArray2D<T>(converted_vals);
}

/** Construct the Atom property from the array of array of QVariants.
     Each array contains the properties for a single CutGroup
    
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::AtomProperty(const QVector< QVector<QVariant> > &values)
                : SireBase::ConcreteProperty<AtomProperty<T>,AtomProp>()
{
    if (values.isEmpty())
        return;
        
    int ngroups = values.count();
    const QVector<QVariant> *values_array = values.constData();
            
    QVector< QVector<T> > tmp_props = QVector< QVector<T> >(ngroups);
    QVector<T> *tmp_props_array = tmp_props.data();
    
    for (int i=0; i<ngroups; ++i)
    {
        const QVector<QVariant> &group_values = values_array[i];
        int nvals = group_values.count();
        
        if (nvals == 0)
            tmp_props_array[i] = QVector<T>();
        else
        {
            QVector<T> converted_vals(nvals);
            T *converted_vals_array = converted_vals.data();
            const QVariant *group_values_array = group_values.constData();

            for (int j=0; j<nvals; ++j)
            {
                const QVariant &value = group_values_array[j];
                AtomProperty<T>::assertCanConvert(value);
                
                converted_vals_array[j] = value.value<T>();
            }
            
            tmp_props_array[i] = converted_vals;
        }
    }
    
    props = PackedArray2D<T>(tmp_props);
}

/** Copy constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::AtomProperty(const AtomProperty<T> &other)
                : SireBase::ConcreteProperty<AtomProperty<T>,AtomProp>(), 
                  props(other.props)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>::~AtomProperty()
{}

/** Copy assignment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>& AtomProperty<T>::operator=(const AtomProperty<T> &other)
{
    props = other.props;
    return *this;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool AtomProperty<T>::operator==(const AtomProperty<T> &other) const
{
    return props == other.props;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool AtomProperty<T>::operator!=(const AtomProperty<T> &other) const
{
    return props != other.props;
}

/** Return the array of properties for the atoms in the CutGroup
    identified by index 'cgidx'
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename PackedArray2D<T>::Array& AtomProperty<T>::operator[](CGIdx cgidx) const
{
    return props.constData()[ cgidx.map(props.count()) ];
}

/** Convert the contained properties into an array of arrays of QVariants.
    There is one array per CutGroup */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QVector< QVector<QVariant> > AtomProperty<T>::toVariant() const
{
    if (props.isEmpty())
        return QVector< QVector<QVariant> >();
        
    int ngroups = props.count();
    const typename PackedArray2D<T>::Array *props_array = props.constData();
    
    QVector< QVector<QVariant> > variant_vals(ngroups);
    variant_vals.squeeze();
    QVector<QVariant> *variant_vals_array = variant_vals.data();
    
    for (int i=0; i<ngroups; ++i)
    {
        const typename PackedArray2D<T>::Array &group_props = props_array[i];
        
        if (group_props.isEmpty())
            variant_vals_array[i] = QVector<QVariant>();
        else
        {
            int nvals = group_props.count();
            const T *group_props_array = group_props.constData();
            
            QVector<QVariant> converted_vals(nvals);
            converted_vals.squeeze();
            QVariant *converted_vals_array = converted_vals.data();
            
            for (int j=0; j<nvals; ++j)
            {
                converted_vals_array[j].setValue<T>( group_props_array[j] );
            }
            
            variant_vals_array[i] = converted_vals;
        }
    }
    
    return variant_vals;
}

/** Assign the values of the properties from the array of QVariants

    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomProperty<T>::assignFrom(const QVector<QVariant> &values) 
{
    this->operator=( AtomProperty<T>(values) );
}

/** Assign the values of the properties from the array of array
    of QVariants
    
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomProperty<T>::assignFrom(const QVector< QVector<QVariant> > &values)
{
    this->operator=( AtomProperty<T>(values) );
}

/** Return the array of properties for the atoms in the CutGroup
    identified by index 'cgidx'
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename PackedArray2D<T>::Array& AtomProperty<T>::at(CGIdx cgidx) const
{
    return this->operator[](cgidx);
}

/** Return the array of properties for the atoms in the CutGroup
    identified by index 'cgidx'
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename PackedArray2D<T>::Array& AtomProperty<T>::get(CGIdx cgidx) const
{
    return this->operator[](cgidx);
}

/** Return the property for the atom at index 'cgatomidx' 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& AtomProperty<T>::operator[](const CGAtomIdx &cgatomidx) const
{
    const typename PackedArray2D<T>::Array &group_props 
                            = this->operator[](cgatomidx.cutGroup());

    return group_props.constData()[ cgatomidx.atom().map(group_props.count()) ];
}

/** Return the property for the atom at index 'cgatomidx' 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& AtomProperty<T>::at(const CGAtomIdx &cgatomidx) const
{
    return this->operator[](cgatomidx);
}

/** Return the property for the atom at index 'cgatomidx' 

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& AtomProperty<T>::get(const CGAtomIdx &cgatomidx) const
{
    return this->operator[](cgatomidx);
}

/** Set the value of the property for the atom at index 'cgatomidx'

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>& AtomProperty<T>::set(const CGAtomIdx &cgatomidx, const T &value)
{
    quint32 cgidx = cgatomidx.cutGroup().map(props.count());
    quint32 atomidx = cgatomidx.atom().map( props.at(cgidx).count() );

    props(cgidx, atomidx) = value;
    
    return *this;
}

/** Set the values for all atoms in the CutGroup at index 'cgidx'

    \throw SireError::incompatible_error
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T>& AtomProperty<T>::set(CGIdx cgidx, const QVector<T> &values)
{
    props.update( cgidx.map(props.count()), values );
    
    return *this;
}

/** Return a const-reference to the PackedArray2D used to store
    all of the atom properties */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const PackedArray2D<T>& AtomProperty<T>::array() const
{
    return props;
}

/** Return a raw pointer to the array of arrays */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename PackedArray2D<T>::Array* AtomProperty<T>::data() const
{
    return props.constData();
}

/** Return a raw pointer to the array of arrays */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const typename PackedArray2D<T>::Array* AtomProperty<T>::constData() const
{
    return props.constData();
}

/** Return a raw pointer to the array of properties for 
    the atoms in the CutGroup at index 'cgidx'
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* AtomProperty<T>::data(CGIdx cgidx) const
{
    return this->at(cgidx).constData();
}

/** Return a raw pointer to the array of properties for 
    the atoms in the CutGroup at index 'cgidx'
    
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T* AtomProperty<T>::constData(CGIdx cgidx) const
{
    return this->at(cgidx).constData();
}

/** Return the number of CutGroups in the molecule */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int AtomProperty<T>::size() const
{
    return props.count();
}

/** Return the number of CutGroups in the molecule */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int AtomProperty<T>::count() const
{
    return props.count();
}

/** Return the number of CutGroups in the molecule */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int AtomProperty<T>::nCutGroups() const
{
    return props.count();
}

/** Return the total number of atoms in the molecule */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int AtomProperty<T>::nAtoms() const
{
    return props.nValues();
}

/** Return the number of atoms in the CutGroup at index 'cgidx'

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int AtomProperty<T>::nAtoms(CGIdx cgidx) const
{
    return this->at(cgidx).count();
}

/** Return whether or not this property is compatible with the
    Molecule whose layout is described in 'molinfo'
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool AtomProperty<T>::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    if (props.nValues() != molinfo.nAtoms())
        return false;

    int ncg = molinfo.nCutGroups();

    if (ncg != props.count())
        return false;
        
    const typename PackedArray2D<T>::Array *props_array = props.constData();
    
    for (CGIdx i(0); i<ncg; ++i)
    {
        if (molinfo.nAtoms(i) != props_array[i].count())
            return false;
    }
    
    return true;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void AtomSelection::assertCompatibleWith(const AtomProperty<T> &prop) const
{
    prop.assertCompatibleWith(this->info());
}

/** Match this property to the passed selection. This returns 
    the property only for the CutGroups that have been selected,
    and with default values for any atoms in those CutGroups that
    have not been selected. This is useful, e.g. for the forcefield
    classes, as this allows an AtomProperty<T> to be returned
    for only the atoms that are selected as part of the forcefield.
    
    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomProperty<T> AtomProperty<T>::matchToSelection(
                                      const AtomSelection &selected_atoms) const
{
    selected_atoms.assertCompatibleWith(*this);

    if (selected_atoms.selectedAll())
        return *this;
    else if (selected_atoms.selectedNone())
        return AtomProperty<T>();
    
    else if (selected_atoms.selectedAllCutGroups())
    {
        PackedArray2D<T> new_props = props;
        
        int ncg = props.count();
        
        for (CGIdx i(0); i<ncg; ++i)
        {
            if (not selected_atoms.selectedAll(i))
            {
                int nats = new_props.at(i).count();
                T *atom_props_array = new_props.data(i);

                const QSet<Index> &atoms = selected_atoms.selectedAtoms(i);
                
                for (Index j(0); j<nats; ++j)
                {
                    if (not atoms.contains(j))
                        atom_props_array[j] = T();
                }
            }
        }
        
        return AtomProperty<T>(new_props);
    }
    else
    {
        QList<CGIdx> cgidxs = selected_atoms.selectedCutGroups();
        
        QVector< QVector<T> > new_props = QVector< QVector<T> >(cgidxs.count());
        QVector<T> *new_props_array = new_props.data();
        
        const typename PackedArray2D<T>::Array *props_array = props.constData();
        
        int n = 0;
        
        foreach (CGIdx i, cgidxs)
        {
            if (selected_atoms.selectedAll(i))
            {
                new_props_array[n] = props_array[i].toQVector();
                ++n;
            }
            else
            {
                const QSet<Index> &atoms = selected_atoms.selectedAtoms(i);
                
                QVector<T> atom_props = props_array[i].toQVector();
                int nats = atom_props.count();
                
                T *atom_props_array = atom_props.data();
                
                for (Index j(0); j<nats; ++j)
                {
                    if (not atoms.contains(j))
                        atom_props_array[j] = T();
                }
                
                new_props_array[n] = atom_props;
                ++n;
            }
        }
        
        return AtomProperty<T>(new_props);
    }
}

}

/** Serialise this property to a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireMol::AtomProperty<T> &prop)
{
    //serialise the base class - this writes the header and version!
    ds << static_cast<const SireMol::AtomProp&>(prop);
    ds << prop.props;
    
    return ds;
}

/** Extract from an binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireMol::AtomProperty<T> &prop)
{
    ds >> static_cast<SireMol::AtomProp&>(prop);
    ds >> prop.props;
        
    return ds;
}

SIRE_END_HEADER

#endif
