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

#ifndef SIREMOL_MOLECULEVIEW_H
#define SIREMOL_MOLECULEVIEW_H

#include <QSharedData>

#include "SireBase/shareddatapointer.hpp"

#include "moleculedata.h"
#include "moleculeinfodata.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class MoleculeView;
}

QDataStream& operator<<(QDataStream&, const SireMol::MoleculeView&);
QDataStream& operator>>(QDataStream&, SireMol::MoleculeView&);

namespace SireMol
{

class MoleculeData;
class AtomSelection;

/** This is the base class of all views of a Molecule. Derived
    classes include Molecule, Segment, Chain, CutGroup, Residue and Atom.
    
    (and the manipulator classes of each of these)

    @author Christopher Woods
*/
class SIREMOL_EXPORT MoleculeView
{

friend QDataStream& ::operator<<(QDataStream&, const MoleculeView&);
friend QDataStream& ::operator>>(QDataStream&, MoleculeView&);

public:
    virtual ~MoleculeView();

    static const char* typeName()
    {
        return "SireMol::MoleculeView";
    }

    /** Return the type name of this view */
    virtual const char* what() const=0;
    
    /** Return a clone of this view */
    virtual MoleculeView* clone() const=0;

    /** Return the MoleculeData that contains all of the information
        about the molecule which this view is showing */
    const MoleculeData& data() const
    {
        return *d;
    }

    /** Return the MoleculeData that contains all of the information
        about the molecule which this view is showing */
    const MoleculeData& constData() const
    {
        return *d;
    }

    /** Return the atoms that are selected as part of this view */
    virtual AtomSelection selection() const=0;

    virtual void update(const MoleculeData &moldata);

    /** Return whether or not this view has the property at key 'key'
         - note that this returns true only if there is a property,
           *and* it fits the view (e.g. is an AtomProperty if this
           is a view of an Atom or group of Atoms) */
    virtual bool hasProperty(const PropertyName &key) const=0;
    
    /** Return whether or not this view has the metadata at metakey 'metakey'
         - note that this returns true only if there is some metadata,
           *and* it fits the view (e.g. is an AtomProperty if this
           is a view of an Atom or group of Atoms) */
    virtual bool hasMetadata(const PropertyName &metakey) const=0;

    /** Return whether or not this view has the metadata at metakey 
        'metakey' for the property at key 'key'
         - note that this returns true only if there is some metadata,
           *and* it fits the view (e.g. is an AtomProperty if this
           is a view of an Atom or group of Atoms) */
    virtual bool hasMetadata(const PropertyName &key, 
                             const PropertyName &metakey) const=0;

    /** Return the list of keys of properties that fit this view,
        e.g. if this is a view of an atom, or group of atoms, then
        this returns the keys of all AtomProperty derived objects */
    virtual QStringList propertyKeys() const=0;
    
    /** Return the list of metakeys of metadata that fit this view,
        e.g. if this is a view of an atom, or group of atoms, then
        this returns the metakeys of all AtomProperty derived objects */
    virtual QStringList metadataKeys() const=0;

    /** Return the list of metakeys of metadata for the 
        property at key 'key' that fit this view,
        e.g. if this is a view of an atom, or group of atoms, then
        this returns the metakeys of all AtomProperty derived objects */
    virtual QStringList metadataKeys(const PropertyName &key) const=0;

    const char* propertyType(const PropertyName &key) const;
    const char* metadataType(const PropertyName &metakey) const;
    
    const char* metadataType(const PropertyName &key,
                             const PropertyName &metakey) const;

    void assertHasProperty(const PropertyName &key) const;
    void assertHasMetadata(const PropertyName &metakey) const;
    void assertHasMetadata(const PropertyName &key,
                           const PropertyName &metakey) const;

    void assertSameMolecule(const MoleculeView &other) const;
    void assertSameMolecule(const MoleculeData &other) const;

    virtual void assertContains(AtomIdx atomidx) const;

protected:
    MoleculeView();
    MoleculeView(const MoleculeData &moldata);
    MoleculeView(const MoleculeView &other);

    MoleculeView& operator=(const MoleculeView &other);

    bool operator==(const MoleculeView &other) const;
    bool operator!=(const MoleculeView &other) const;

    template<class Index, class PropType, class T>
    static void setProperty(MoleculeData &data,
                            const QString &key, 
                            const Index &idx, const T &value);

    template<class Index, class PropType, class T>
    static void setMetadata(MoleculeData &data,
                            const QString &metakey, 
                            const Index &idx, const T &value);

    template<class Index, class PropType, class T>
    static void setMetadata(MoleculeData &data,
                            const QString &key, const QString &metakey,
                            const Index &idx, const T &value);
                            

    /** Shared pointer to the raw data of the molecule */
    SireBase::SharedDataPointer<MoleculeData> d;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

template<class Index, class PropType, class T>
SIRE_OUTOFLINE_TEMPLATE
void MoleculeView::setProperty(MoleculeData &data,
                               const QString &key, 
                               const Index &idx, const T &value)
{
    PropType props;
    
    if (data.hasProperty(key))
    {
        //take the property to prevent unnecessary copying caused
        //by implicit sharing of the property
        Property old_property = data.takeProperty(key);
        props = old_property->asA<PropType>();
    }
    else
        props = PropType(data.info());
        
    props.set(idx, value);
    
    data.setProperty(key, props);
}

template<class Index, class PropType, class T>
SIRE_OUTOFLINE_TEMPLATE
void MoleculeView::setMetadata(MoleculeData &data,
                               const QString &metakey, 
                               const Index &idx, const T &value)
{
    PropType props;
    
    if (data.hasMetadata(metakey))
    {
        //take the metadata to prevent unnecessary copying caused
        //by implicit sharing of the property
        Property old_metadata = data.takeMetadata(metakey);
        props = old_metadata->asA<PropType>();
    }
    else
        props = PropType(data);
        
    props.set(idx, value);
    
    data.setMetadata(metakey, props);
}

template<class Index, class PropType, class T>
SIRE_OUTOFLINE_TEMPLATE
void MoleculeView::setMetadata(MoleculeData &data,
                               const QString &key, const QString &metakey,
                               const Index &idx, const T &value)
{
    PropType props;
    
    if (data.hasMetadata(key, metakey))
    {
        //take the metadata to prevent unnecessary copying caused
        //by implicit sharing of the property
        Property old_metadata = data.takeMetadata(key, metakey);
        props = old_metadata->asA<PropType>();
    }
    else
        props = PropType(data);
        
    props.set(idx, value);
    
    data.setMetadata(key, metakey, props);
}

namespace detail
{

template<class Prop, class Idx, class V>
QList<V> get_property(const MoleculeData &moldata,
                      const QList<Idx> &idxs,
                      const PropertyName &key)
{
    QList<V> props;
        
    const Property &property = moldata.property(key);
    
    const Prop &prop = property->asA<Prop>();
        
    foreach (Idx idx, idxs)
    {
        props.append( prop.at(idx) );
    }
        
    return props;
}
    
template<class Prop, class Idx, class V>
QList<V> get_metadata(const MoleculeData &moldata,
                      const QList<Idx> &idxs,
                      const PropertyName &metakey)
{
    QList<V> props;
    
    const Property &property = moldata.metadata(metakey);
    const Prop &prop = property->asA<Prop>();
    
    foreach (Idx idx, idxs)
    {
        props.append( prop.at(idx) );
    }
    
    return props;
}

template<class Prop, class Idx, class V>
QList<V> get_metadata(const MoleculeData &moldata,
                      const QList<Idx> &idxs,
                      const PropertyName &key,
                      const PropertyName &metakey)
{
    QList<V> props;
    
    const Property &property = moldata.metadata(key,metakey);
    const Prop &prop = property->asA<Prop>();
    
    foreach (Idx idx, idxs)
    {
        props.append( prop.at(idx) );
    }
    
    return props;
}

template<class Prop, class Idx, class V>
void set_property(MoleculeData &moldata,
                  const QList<Idx> &idxs,
                  const QString &key,
                  const QList<V> &values)
{
    Prop props;
    
    if (moldata.hasProperty(key))
        props = moldata.property(key);
    else
        props = Prop(moldata.info());
        
    for (int i=0; i<idxs.count(); ++i)
    {
        props.set(idxs[i], values[i]);
    }
    
    moldata.setProperty(key, props);
}

template<class Prop, class V, class Idx>
void set_metadata(MoleculeData &moldata,
                  const QList<Idx> &idxs,
                  const QString &metakey,
                  const QList<V> &values)
{
    Prop props;
    
    if (moldata.hasMetadata(metakey))
        props = moldata.metadata(metakey);
    else
        props = Prop(moldata.info());
        
    for (int i=0; i<idxs.count(); ++i)
    {
        props.set(idxs[i], values[i]);
    }
    
    moldata.setMetadata(metakey, props);
}

template<class Prop, class Idx, class V>
void set_metadata(MoleculeData &moldata,
                  const QList<Idx> &idxs,
                  const QString &key, const QString &metakey,
                  const QList<V> &values)
{
    Prop props;
    
    if (moldata.hasMetadata(key,metakey))
        props = moldata.metadata(key,metakey);
    else
        props = Prop(moldata.info());
        
    for (int i=0; i<idxs.count(); ++i)
    {
        props.set(idxs[i], values[i]);
    }
    
    moldata.setMetadata(key, metakey, props);
}

template<class Prop, class Idx, class V>
void set_property(MoleculeData &moldata,
                  const QList<Idx> &idxs,
                  const QString &key,
                  const V &value)
{
    Prop props;
    
    if (moldata.hasProperty(key))
        props = moldata.property(key);
    else
        props = Prop(moldata.info());
        
    foreach (Idx idx, idxs)
    {
        props.set(idx, value);
    }
    
    moldata.setProperty(key, props);
}

template<class Prop, class Idx,class V>
void set_metadata(MoleculeData &moldata,
                  const QList<Idx> &idxs,
                  const QString &metakey,
                  const V &value)
{
    Prop props;
    
    if (moldata.hasMetadata(metakey))
        props = moldata.metadata(metakey);
    else
        props = Prop(moldata.info());
        
    foreach (Idx idx, idxs)
    {
        props.set(idx, value);
    }
    
    moldata.setMetadata(metakey, props);
}

template<class Prop, class Idx, class V>
void set_metadata(MoleculeData &moldata,
                  const QList<Idx> &idxs,
                  const QString &key, const QString &metakey,
                  const V &value)
{
    Prop props;
    
    if (moldata.hasMetadata(key,metakey))
        props = moldata.metadata(key,metakey);
    else
        props = Prop(moldata.info());
        
    foreach (Idx idx, idxs)
    {
        props.set(idx, value);
    }
    
    moldata.setMetadata(key, metakey, props);
}

}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_EXPOSE_CLASS( SireMol::MoleculeView )

SIRE_END_HEADER

#endif
