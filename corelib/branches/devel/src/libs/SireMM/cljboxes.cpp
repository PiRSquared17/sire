/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#include "SireMM/cljboxes.h"

#include "SireVol/aabox.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"

#include <QElapsedTimer>
#include <QDebug>

using namespace SireMM;
using namespace SireMaths;
using namespace SireVol;
using namespace SireStream;

///////////
/////////// Implementation of CLJBox
///////////

static const RegisterMetaType<CLJBox> r_cljbox(NO_ROOT);

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const CLJBox &box)
{
    writeHeader(ds, r_cljbox, 1);
    
    SharedDataStream sds(ds);
    sds << box.atms;
    
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, CLJBox &box)
{
    VersionID v = readHeader(ds, r_cljbox);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> box.atms;
    }
    else
        throw version_error(v, "1", r_cljbox, CODELOC);
    
    return ds;
}

/** Null constructor */
CLJBox::CLJBox() : QSharedData()
{}

/** Copy constructor */
CLJBox::CLJBox(const CLJBox &other) : QSharedData(), atms(other.atms)
{}

/** Destructor */
CLJBox::~CLJBox()
{}

/** Copy assignment operator */
CLJBox& CLJBox::operator=(const CLJBox &other)
{
    atms = other.atms;
    return *this;
}

/** Comparison operator */
bool CLJBox::operator==(const CLJBox &other) const
{
    return atms == other.atms;
}

/** Comparison operator */
bool CLJBox::operator!=(const CLJBox &other) const
{
    return not operator==(other);
}

const char* CLJBox::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJBox>() );
}

const char* CLJBox::what() const
{
    return CLJBox::typeName();
}

///////////
/////////// Implementation of CLJBoxPtr
///////////

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const CLJBoxPtr &ptr)
{
    SharedDataStream sds(ds);
    sds << ptr.d;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, CLJBoxPtr &ptr)
{
    SharedDataStream sds(ds);
    sds >> ptr.d;
    return ds;
}

static const QSharedDataPointer<CLJBox> *null_box = 0;

/** Constructor */
CLJBoxPtr::CLJBoxPtr()
{
    if (null_box == 0)
    {
        null_box = new QSharedDataPointer<CLJBox>(new CLJBox());
    }
    
    d = *null_box;
}

/** Construct from the passed box. This takes over ownership of the pointer,
    and will delete it once there are no more CLJBoxPtrs pointing to this box */
CLJBoxPtr::CLJBoxPtr(CLJBox *box)
{
    if (box == 0)
    {
        if (null_box == 0)
        {
            null_box = new QSharedDataPointer<CLJBox>(new CLJBox());
        }
        
        d = *null_box;
    }
    else
    {
        d = QSharedDataPointer<CLJBox>(box);
    }
}

/** Copy constructor */
CLJBoxPtr::CLJBoxPtr(const CLJBoxPtr &other) : d(other.d)
{}

/** Destructor */
CLJBoxPtr::~CLJBoxPtr()
{
    //this will delete the box if it is the last
    //reference to the box
}

/** Copy assignment operator */
CLJBoxPtr& CLJBoxPtr::operator=(const CLJBoxPtr &other)
{
    d = other.d;
    return *this;
}

/** Comparison operator */
bool CLJBoxPtr::operator==(const CLJBoxPtr &other) const
{
    return d.constData() == other.d.constData() or
           *(d.constData()) == *(other.d.constData());
}

/** Comparison operator */
bool CLJBoxPtr::operator!=(const CLJBoxPtr &other) const
{
    return not operator==(other);
}

///////////
/////////// Implementation of CLJBoxIndex
///////////

static const RegisterMetaType<CLJBoxIndex> r_cljboxindex(NO_ROOT);

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const CLJBoxIndex &index)
{
    writeHeader(ds, r_cljboxindex, 1);
    ds << index.v.index.ii << index.v.index.jj << index.v.index.kk << index.v.index.idx;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, CLJBoxIndex &index)
{
    VersionID v = readHeader(ds, r_cljboxindex);
    
    if (v == 1)
    {
        ds >> index.v.index.ii >> index.v.index.jj
           >> index.v.index.kk >> index.v.index.idx;
    }
    else
        throw version_error(v, "1", r_cljboxindex, CODELOC);
    
    return ds;
}

/** Null constructor */
CLJBoxIndex::CLJBoxIndex()
{
    v.val = 0;
}

/** Construct the index of the box at index i,j,k with (optionally supplied)
    index of a particular atom in the box */
CLJBoxIndex::CLJBoxIndex(qint16 ii, qint16 jj, qint16 kk, quint16 atom_idx)
{
    v.index.ii = ii;
    v.index.jj = jj;
    v.index.kk = kk;
    v.index.idx = atom_idx;
}

/** Copy constructor */
CLJBoxIndex::CLJBoxIndex(const CLJBoxIndex &other)
{
    v.val = other.v.val;
}

/** Destructor */
CLJBoxIndex::~CLJBoxIndex()
{}

/** Copy assignment operator */
CLJBoxIndex& CLJBoxIndex::operator=(const CLJBoxIndex &other)
{
    v.val = other.v.val;
    return *this;
}

const char* CLJBoxIndex::typeName()
{
    return "SireMM::CLJBoxIndex";
}

const char* CLJBoxIndex::what() const
{
    return CLJBoxIndex::typeName();
}

/** Return the AABox that describes this box (for a given box length of 'box_length') */
AABox CLJBoxIndex::box(Length box_length) const
{
    Vector length(box_length.value());
    
    Vector origin( i() * length.x(),
                   j() * length.y(),
                   k() * length.z() );
    
    return AABox(origin, 0.5*length);
}

/** Return a copy of this index that contains only the box index (not the atom index) */
CLJBoxIndex CLJBoxIndex::boxOnly() const
{
    return CLJBoxIndex( i(), j(), k() );
}

/** Create the index for the box that contains the point 'x,y,z' in a set of boxes
    of length 1 / inv_box_length */
CLJBoxIndex CLJBoxIndex::createWithInverseBoxLength(float x, float y, float z, float inv_length)
{
    int i = (x * inv_length) + 0.5;
    int j = (y * inv_length) + 0.5;
    int k = (z * inv_length) + 0.5;
    
    const int min16 = std::numeric_limits<qint16>::min();
    const int max16 = std::numeric_limits<qint16>::max();
    
    if (i < min16 or i > max16 or
        j < min16 or j > max16 or
        k < min16 or k > max16)
    {
        throw SireError::program_bug( QObject::tr(
                "It is not possible to get the index of the box containing the point "
                "(%1,%2,%3) when the boxes have length %4 A, as the index lies outside "
                "the bounds of a 16bit integer (%5 to %6), (%7,%8,%9).")
                    .arg(x).arg(y).arg(z)
                    .arg(1 / inv_length)
                    .arg(min16).arg(max16)
                    .arg(i).arg(j).arg(k), CODELOC );
    }
    
    return CLJBoxIndex(i, j, k);
}

/** Create the index for the box that contains the point 'x,y,z' in a set of boxes
    of length 1 / inv_box_length */
CLJBoxIndex CLJBoxIndex::createWithInverseBoxLength(const Vector &coords, float inv_length)
{
    return createWithInverseBoxLength(coords.x(), coords.y(), coords.z(), inv_length);
}

/** Create the index for the box that contains the point 'x,y,z' in a set of boxes
    of length 'box_length' */
CLJBoxIndex CLJBoxIndex::createWithBoxLength(float x, float y, float z, Length box_length)
{
    return CLJBoxIndex::createWithInverseBoxLength(x, y, z, 1.0/box_length.value());
}

/** Create the index for the box that contains the point 'x,y,z' in a set of boxes
    of length 'box_length' */
CLJBoxIndex CLJBoxIndex::createWithBoxLength(const Vector &coords, Length box_length)
{
    return CLJBoxIndex::createWithInverseBoxLength(coords.x(), coords.y(),
                                                   coords.z(), 1.0/box_length.value());
}

///////////
/////////// Implementation of CLJBoxes
///////////

static const RegisterMetaType<CLJBoxes> r_cljboxes(NO_ROOT);

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const CLJBoxes &boxes)
{
    writeHeader(ds, r_cljboxes, 1);
    
    SharedDataStream sds(ds);
    
    sds << boxes.bxs << boxes.box_length;
    
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, CLJBoxes &boxes)
{
    VersionID v = readHeader(ds, r_cljboxes);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> boxes.bxs >> boxes.box_length;
    }
    else
        throw version_error(v, "1", r_cljboxes, CODELOC);
    
    return ds;
}

const double default_box_length = 5;

/** Null constructor */
CLJBoxes::CLJBoxes() : box_length(default_box_length)
{}

/** Construct from the passed set of atoms */
void CLJBoxes::constructFrom(const CLJAtoms &atoms)
{
    QElapsedTimer t;
    t.start();
    
    const QVector<MultiFloat> &xm = atoms.x();
    const QVector<MultiFloat> &ym = atoms.y();
    const QVector<MultiFloat> &zm = atoms.z();
    const QVector<MultiFloat> &qm = atoms.q();
    const QVector<MultiFloat> &sigm = atoms.sigma();
    const QVector<MultiFloat> &epsm = atoms.epsilon();
    const QVector<MultiInt> &idm = atoms.ID();
    
    const MultiFloat *xa = xm.constData();
    const MultiFloat *ya = ym.constData();
    const MultiFloat *za = zm.constData();
    const MultiFloat *qa = qm.constData();
    const MultiFloat *siga = sigm.constData();
    const MultiFloat *epsa = epsm.constData();
    const MultiInt *ida = idm.constData();
    
    for (int i=0; i<xm.count(); ++i)
    {
        for (int ii=0; ii<MultiFloat::count(); ++ii)
        {
            if (ida[i][ii] != 0)
            {
                //box up this atom
                const float x = xa[i][ii];
                const float y = ya[i][ii];
                const float z = za[i][ii];
            }
        }
    }
    
    quint64 ns = t.elapsed();
    qDebug() << "Boxing up" << atoms.count() << "atoms took"
             << (0.000001*ns) << "ms";
}

/** Box up the passed set of atoms into boxes of the default box size
    (5 angstroms) */
CLJBoxes::CLJBoxes(const CLJAtoms &atoms) : box_length(default_box_length)
{
    constructFrom(atoms);
}

/** Box up the passed set of atoms into boxes of specified box size */
CLJBoxes::CLJBoxes(const CLJAtoms &atoms, Length box_size)
         : box_length(box_size.value())
{
    if (box_length <= 0)
        box_length = default_box_length;
    
    constructFrom(atoms);
}

/** Copy constructor */
CLJBoxes::CLJBoxes(const CLJBoxes &other)
         : bxs(other.bxs), box_length(other.box_length)
{}

/** Destructor */
CLJBoxes::~CLJBoxes()
{}

/** Copy assignment operator */
CLJBoxes& CLJBoxes::operator=(const CLJBoxes &other)
{
    if (this != &other)
    {
        bxs = other.bxs;
        box_length = other.box_length;
    }
    
    return *this;
}

/** Comparison operator */
bool CLJBoxes::operator==(const CLJBoxes &other) const
{
    return box_length == other.box_length and bxs == other.bxs;
}

/** Comparison operator */
bool CLJBoxes::operator!=(const CLJBoxes &other) const
{
    return not operator==(other);
}

/** Return the indicies of all occupied boxes */
QVector<CLJBoxIndex> CLJBoxes::occupiedBoxIndicies() const
{
    return bxs.keys().toVector();
}

/** Return the box at index 'index' */
CLJBox CLJBoxes::boxAt(const CLJBoxIndex &index) const
{
    return bxs.value( index.boxOnly() ).read();
}

/** Return the AABox that describes the boundary of the box at index 'index' */
AABox CLJBoxes::boxDimensionsAt(const CLJBoxIndex &index) const
{
    return index.box( Length(box_length) );
}

/** Return the box that contains the point with coordinates 'coords' */
CLJBox CLJBoxes::boxAt(const Vector &coords) const
{
    return this->boxAt( CLJBoxIndex::createWithInverseBoxLength(coords, 1.0/box_length) );
}

/** Return the dimensions of the box that contains the point with coordinates 'coords' */
AABox CLJBoxes::boxDimensionsAt(const Vector &coords) const
{
    return this->boxDimensionsAt( CLJBoxIndex::createWithInverseBoxLength(coords, 1.0/box_length));
}

/** Return an array containing all occupied boxes. This is in the same order
    as the box dimensions returned by 'boxDimensions' */
QVector<CLJBox> CLJBoxes::boxes() const
{
    if (bxs.isEmpty())
        return QVector<CLJBox>();

    QVector<CLJBox> b( bxs.count() );
    CLJBox *ba = b.data();
    
    int idx = 0;
    
    for (QMap<CLJBoxIndex,CLJBoxPtr>::const_iterator it = bxs.constBegin();
         it != bxs.constEnd();
         ++it)
    {
        ba[idx] = it.value().read();
        idx += 1;
    }
    
    return b;
}

/** Return an array containing all occupied box dimensions. This is in the same order
    as the boxes returned by 'boxes' */
QVector<AABox> CLJBoxes::boxDimensions() const
{
    if (bxs.isEmpty())
        return QVector<AABox>();

    QVector<AABox> b( bxs.count() );
    AABox *ba = b.data();
    
    int idx = 0;
    
    for (QMap<CLJBoxIndex,CLJBoxPtr>::const_iterator it = bxs.constBegin();
         it != bxs.constEnd();
         ++it)
    {
        ba[idx] = it.key().box(Length(box_length));
        idx += 1;
    }
    
    return b;
}

/** Return all of the atoms in all of the boxes (these may
    be returned with a lot of padding) */
CLJAtoms CLJBoxes::atoms() const
{
    CLJAtoms atms;
    
    for (QMap<CLJBoxIndex,CLJBoxPtr>::const_iterator it = bxs.constBegin();
         it != bxs.constEnd();
         ++it)
    {
        atms += it.value().read().atoms();
    }
    
    return atms;
}

/** Return the length of each side of each box */
Length CLJBoxes::length() const
{
    return Length(box_length);
}
