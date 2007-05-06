/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREMOVE_UNIFORMSAMPLER_H
#define SIREMOVE_UNIFORMSAMPLER_H

#include "sampler.h"

#include "SireMol/moleculeid.h"
#include "SireMol/moleculegroupid.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class UniformSampler;
}

QDataStream& operator<<(QDataStream&, const SireMove::UniformSampler&);
QDataStream& operator>>(QDataStream&, SireMove::UniformSampler&);

namespace SireMove
{

using SireMol::MoleculeID;
using SireMol::MoleculeGroupID;

/** This class is used to pick contained Molecules, Residues
    or Atoms at random (uniformly distributed).

    @author Christopher Woods
*/
class SIREMOVE_EXPORT UniformSampler : public SamplerBase
{

friend QDataStream& ::operator<<(QDataStream&, const UniformSampler&);
friend QDataStream& ::operator>>(QDataStream&, UniformSampler&);

public:
    UniformSampler();
    UniformSampler(const RanGenerator &rangenerator);

    UniformSampler(const UniformSampler &other);

    ~UniformSampler();

    UniformSampler& operator=(const UniformSampler &other);

    static const char* typeName()
    {
        return "SireMove::UniformSampler";
    }

    const char* what() const
    {
        return UniformSampler::typeName();
    }

    UniformSampler* clone() const
    {
        return new UniformSampler(*this);
    }

    tuple<PartialMolecule,double> sample(const MoleculeGroup &group);

    double probabilityOf(const PartialMolecule &molecule,
                         const MoleculeGroup &group);

protected:
    bool _pvt_isEqual(const PropertyBase &other) const;

private:
    void updateFrom(const MoleculeGroup &group);

    /** List of IDs of all of the molecules that could be selected */
    QList<MoleculeID> molids;
    
    /** The ID of the MoleculeGroup from which the last molecule
        was chosen */
    MoleculeGroupID groupid;
    
    /** The major version number of the group from which the last molecule
        was chosen */
    quint32 majver;
};

}

Q_DECLARE_METATYPE(SireMove::UniformSampler);

SIRE_END_HEADER

#endif
