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

#ifndef SIREMOVE_RIGIDBODYMC_H
#define SIREMOVE_RIGIDBODYMC_H

#include "montecarlo.h"
#include "sampler.h"
#include "getpoint.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class RigidBodyMC;
}

QDataStream& operator<<(QDataStream&, const SireMove::RigidBodyMC&);
QDataStream& operator>>(QDataStream&, SireMove::RigidBodyMC&);

namespace SireMol
{
class MoleculeGroup;
class PartialMolecule;
}

namespace SireMove
{

class Sampler;

using SireMol::MoleculeGroup;
using SireMol::PartialMolecule;

/** This class implements a rigid body Monte Carlo move that
    may be applied to a random molecule within a MoleculeGroup

    @author Christopher Woods
*/
class SIREMOVE_EXPORT RigidBodyMC 
        : public SireBase::ConcreteProperty<RigidBodyMC,MonteCarlo>
{

friend QDataStream& ::operator<<(QDataStream&, const RigidBodyMC&);
friend QDataStream& ::operator>>(QDataStream&, RigidBodyMC&);

public:
    RigidBodyMC(const PropertyMap &map = PropertyMap());

    RigidBodyMC(const MoleculeGroup &molgroup,
                const PropertyMap &map = PropertyMap());
    RigidBodyMC(const Sampler &sampler,
                const PropertyMap &map = PropertyMap());

    RigidBodyMC(const RigidBodyMC &other);

    ~RigidBodyMC();

    RigidBodyMC& operator=(const RigidBodyMC &other);

    static const char* typeName();

    bool operator==(const RigidBodyMC &other) const;
    bool operator!=(const RigidBodyMC &other) const;

    QString toString() const;

    void setSampler(const Sampler &sampler);
    void setSampler(const MoleculeGroup &molgroup);

    const Sampler& sampler() const;
    const MoleculeGroup& moleculeGroup() const;

    void setGenerator(const RanGenerator &rangenerator);

    void setMaximumTranslation(SireUnits::Dimension::Length max_translation);
    void setMaximumRotation(SireUnits::Dimension::Angle max_rotation);

    void setCenterOfRotation(const GetPoint &center_function);

    void setReflectionSphere(Vector sphere_center,
                             SireUnits::Dimension::Length sphere_radius);

    bool usesReflectionMoves() const;
    Vector reflectionSphereCenter() const;
    SireUnits::Dimension::Length reflectionSphereRadius() const;

    void setSynchronisedTranslation(bool on);
    void setSynchronisedRotation(bool on);
    void setSharedRotationCenter(bool on);
    
    bool synchronisedTranslation() const;
    bool synchronisedRotation() const;
    bool sharedRotationCenter() const;

    SireUnits::Dimension::Length maximumTranslation() const;
    SireUnits::Dimension::Angle maximumRotation() const;

    const GetPoint& centerOfRotation() const;

    void move(System &system, int nmoves, bool record_stats=true);

protected:
    void _pvt_setTemperature(const SireUnits::Dimension::Temperature &temperature);

private:
    void performMove(System &system,
                     double &old_bias, double &new_bias,
                     const SireBase::PropertyMap &map);

    /** The sampler used to select random molecules for the move */
    SamplerPtr smplr;

    /** The function used to get the center of rotation of each molecule */
    GetPointPtr center_function;

    /** The maximum translation */
    double adel;

    #ifndef SKIP_BROKEN_GCCXML_PARTS
    /** The maximum rotation */
    SireUnits::Dimension::Angle rdel;
    #endif
    
    /** The center of the reflection sphere */
    SireMaths::Vector reflect_center;

    /** The radius of the reflection sphere */
    double reflect_radius;

    /** Whether or not to reflect the moves */
    bool reflect_moves;

    /** Whether or not to synchronise the translation of all views */
    bool sync_trans;
    
    /** Whether or not to synchronise the rotation of all views */
    bool sync_rot;

    /** Whether or not to use the same center of rotation for
        all views - this only applies when synchronised rotation
        is on */
    bool common_center;
};

}

Q_DECLARE_METATYPE(SireMove::RigidBodyMC);

SIRE_EXPOSE_CLASS( SireMove::RigidBodyMC )

SIRE_END_HEADER

#endif
