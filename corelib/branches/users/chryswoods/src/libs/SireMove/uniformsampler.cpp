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

#include <QMutex>

#include "uniformsampler.h"

#include "SireSystem/system.h"

#include "SireMol/partialmolecule.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

using namespace SireMol;
using namespace SireMove;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<UniformSampler> r_sampler;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds,
                                        const UniformSampler &sampler)
{
    writeHeader(ds, r_sampler, 1)
          << static_cast<const Sampler&>(sampler);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, UniformSampler &sampler)
{
    VersionID v = readHeader(ds, r_sampler);

    if (v == 1)
    {
        ds >> static_cast<Sampler&>(sampler);
    }
    else
        throw version_error(v, "1", r_sampler, CODELOC);

    return ds;
}

/** Constructor */
UniformSampler::UniformSampler() 
               : ConcreteProperty<UniformSampler,Sampler>()
{}

static SharedPolyPointer<UniformSampler> shared_null;

const UniformSampler& Sampler::null()
{
    if (shared_null.constData() == 0)
    {
        QMutexLocker lkr( SireBase::globalLock() );
        
        if (shared_null.constData() == 0)
            shared_null = new UniformSampler();
    }
    
    return *(shared_null.constData());
}

/** Constructor a sampler that chooses views at random from the 
    passed molecule group */
UniformSampler::UniformSampler(const MoleculeGroup &molgroup)
               : ConcreteProperty<UniformSampler,Sampler>(molgroup)
{}

/** Copy constructor */
UniformSampler::UniformSampler(const UniformSampler &other)
               : ConcreteProperty<UniformSampler,Sampler>(other)
{}

/** Destructor */
UniformSampler::~UniformSampler()
{}

/** Copy assignment */
UniformSampler& UniformSampler::operator=(const UniformSampler &other)
{
    Sampler::operator=(other);

    return *this;
}

/** Return a random view molecule from the molecule group, together with
    the probability of choosing that view. */
tuple<PartialMolecule,double> UniformSampler::sample() const
{
    //how many molecule views are there in this molecule group?
    int nviews = group().nViews();

    if (nviews == 0)
        throw SireMol::missing_molecule( QObject::tr(
            "The MoleculeGroup is empty, so we can't choose a molecule!"),
                CODELOC );
    
    else if (nviews == 1)
        return tuple<PartialMolecule,double>(group().viewAt(0), 1);

    //choose a random view
    quint32 i = generator().randInt(nviews-1);
    
    //return the matching molecule
    return tuple<PartialMolecule,double>(group().viewAt(i), 1.0 / nviews);
}

/** Return the probability of selecting the view in 'molview' from
    the system 'system'. A probability of zero is returned if
    this view cannot be chosen from the molecule group. */
double UniformSampler::probabilityOf(const PartialMolecule &molecule) const
{
    if (group().contains(molecule))
        return 1.0 / group().nViews();
        
    else
        return 0;
}

const char* UniformSampler::typeName()
{
    return QMetaType::typeName( qMetaTypeId<UniformSampler>() );
}

UniformSampler* UniformSampler::clone() const
{
    return new UniformSampler(*this);
}

