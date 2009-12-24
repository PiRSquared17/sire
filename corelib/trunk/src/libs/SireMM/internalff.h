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

#ifndef SIREMM_INTERNALFF_H
#define SIREMM_INTERNALFF_H

#include "SireBase/propertymap.h"

#include "SireFF/g1ff.h"
#include "SireFF/ff3d.h"
#include "SireFF/forcetable.h"

#include "SireFF/detail/ffmolecules3d.h"

#include "internalcomponent.h"
#include "internalparameters.h"

namespace SireMM
{
class InternalFF;
}

QDataStream &operator<<(QDataStream&, const SireMM::InternalFF&);
QDataStream &operator>>(QDataStream&, SireMM::InternalFF&);

namespace SireVol
{
class CoordGroup;
}

namespace SireMM
{

using SireFF::ForceField;
using SireFF::G1FF;
using SireFF::FF3D;
using SireFF::ForceTable;
using SireFF::MolForceTable;
using SireFF::FFComponent;

using SireMol::MolNum;
using SireMol::PartialMolecule;
using SireMol::MoleculeGroup;

using SireVol::CoordGroup;

using SireBase::Property;
using SireBase::Properties;
using SireBase::PropertyMap;
using SireBase::PropertyName;

/** This class provides the default name of the 
    property that contains the bond parameters */
class SIREMM_EXPORT BondParameterName
{
public:
    BondParameterName()
    {}
    
    ~BondParameterName()
    {}
    
    const PropertyName& bond() const
    {
        return bond_param;
    }
    
private:
    static PropertyName bond_param;
};

/** This class provides the default name of the 
    property that contains the angle parameters */
class SIREMM_EXPORT AngleParameterName
{
public:
    AngleParameterName()
    {}
    
    ~AngleParameterName()
    {}
    
    const PropertyName& angle() const
    {
        return angle_param;
    }
    
private:
    static PropertyName angle_param;
};

/** This class provides the default name of the 
    property that contains the dihedral parameters */
class SIREMM_EXPORT DihedralParameterName
{
public:
    DihedralParameterName()
    {}
    
    ~DihedralParameterName()
    {}
    
    const PropertyName& dihedral() const
    {
        return dihedral_param;
    }
    
private:
    static PropertyName dihedral_param;
};

/** This class provides the default name of the 
    property that contains the improper parameters */
class SIREMM_EXPORT ImproperParameterName
{
public:
    ImproperParameterName()
    {}
    
    ~ImproperParameterName()
    {}
    
    const PropertyName& improper() const
    {
        return improper_param;
    }
    
private:
    static PropertyName improper_param;
};

/** This class provides the default name of the 
    property that contains the Urey-Bradley parameters */
class SIREMM_EXPORT UreyBradleyParameterName
{
public:
    UreyBradleyParameterName()
    {}
    
    ~UreyBradleyParameterName()
    {}
    
    const PropertyName& ureyBradley() const
    {
        return ub_param;
    }
    
private:
    static PropertyName ub_param;
};

/** This class provides the default name of the 
    property that contains the stretch-stretch parameters */
class SIREMM_EXPORT StretchStretchParameterName
{
public:
    StretchStretchParameterName()
    {}
    
    ~StretchStretchParameterName()
    {}
    
    const PropertyName& stretchStretch() const
    {
        return ss_param;
    }
    
private:
    static PropertyName ss_param;
};

/** This class provides the default name of the 
    property that contains the stretch-bend parameters */
class SIREMM_EXPORT StretchBendParameterName
{
public:
    StretchBendParameterName()
    {}
    
    ~StretchBendParameterName()
    {}
    
    const PropertyName& stretchBend() const
    {
        return sb_param;
    }
    
private:
    static PropertyName sb_param;
};

/** This class provides the default name of the 
    property that contains the bend-bend parameters */
class SIREMM_EXPORT BendBendParameterName
{
public:
    BendBendParameterName()
    {}
    
    ~BendBendParameterName()
    {}
    
    const PropertyName& bendBend() const
    {
        return bb_param;
    }
    
private:
    static PropertyName bb_param;
};

/** This class provides the default name of the 
    property that contains the stretch-bend-torsion parameters */
class SIREMM_EXPORT StretchBendTorsionParameterName
{
public:
    StretchBendTorsionParameterName()
    {}
    
    ~StretchBendTorsionParameterName()
    {}
    
    const PropertyName& stretchBendTorsion() const
    {
        return sbt_param;
    }
    
private:
    static PropertyName sbt_param;
};

/** This class provides the default name of the properties
    that contain the bond, angle, dihedral and Urey-Bradley parameters */
class SIREMM_EXPORT InternalParameterNames 
            : public BondParameterName, public AngleParameterName,
              public DihedralParameterName, 
              public ImproperParameterName, public UreyBradleyParameterName,
              public StretchStretchParameterName,
              public StretchBendParameterName,
              public BendBendParameterName,
              public StretchBendTorsionParameterName
{
public:
    InternalParameterNames() 
            : BondParameterName(), AngleParameterName(),
              DihedralParameterName(), 
              ImproperParameterName(), UreyBradleyParameterName(),
              StretchStretchParameterName(), StretchBendParameterName(),
              BendBendParameterName(), StretchBendTorsionParameterName()
    {}
    
    ~InternalParameterNames()
    {}
};

/** This class provides the default name of the properties
    that contain the internal and 3D coordinates properties */
class SIREMM_EXPORT InternalParameterNames3D 
                    : public InternalParameterNames,
                      public SireFF::detail::Coords3DParameterName
{
public:
    InternalParameterNames3D() : InternalParameterNames(),
                                 SireFF::detail::Coords3DParameterName()
    {}
    
    ~InternalParameterNames3D()
    {}
};

/** This potential provides energies and forces caused by 
    molecular mechanics style internal intramolecular 
    terms, e.g. bond, angle, dihedral, urey bradley
    
    @author Christopher Woods
*/
class SIREMM_EXPORT InternalPotential
{
public:
    typedef InternalEnergy Energy;
    typedef Energy::Components Components;
    
    typedef InternalParameterNames3D ParameterNames;
    
    typedef InternalParameters3D Parameters;

    typedef SireFF::detail::FFMolecule<InternalPotential> Molecule;
    typedef SireFF::detail::FFMolecules<InternalPotential> Molecules;
    typedef Molecules::ChangedMolecule ChangedMolecule;
    
    ~InternalPotential();
    
    static const InternalSymbols& symbols();
    
    static const ParameterNames& parameters();

    InternalPotential::Parameters
    getParameters(const PartialMolecule &molecule,
                  const PropertyMap &map) const;

    InternalPotential::Parameters
    updateParameters(const InternalPotential::Parameters &old_params,
                     const PartialMolecule &old_molecule,
                     const PartialMolecule &new_molecule,
                     const PropertyMap &map) const;

    InternalPotential::Parameters
    updateParameters(const InternalPotential::Parameters &old_params,
                     const PartialMolecule &old_molecule,
                     const PartialMolecule &new_molecule,
                     const PropertyMap &old_map,
                     const PropertyMap &new_map) const;

    InternalPotential::Molecule
    parameterise(const PartialMolecule &molecule,
                 const PropertyMap &map) const;
                 
    InternalPotential::Molecules
    parameterise(const MoleculeGroup &molecules,
                 const PropertyMap &map) const;
    
protected:
    InternalPotential(bool isstrict=false);
    InternalPotential(const InternalPotential &other);
    
    InternalPotential& operator=(const InternalPotential &other);

    void calculateEnergy(const InternalPotential::Molecule &molecule,
                         InternalPotential::Energy &energy,
                         double scale_energy=1) const;
    
    void calculateForce(const InternalPotential::Molecule &molecule,
                        MolForceTable &forces,
                        double scale_force=1) const;

    void calculateBondForce(const InternalPotential::Molecule &molecule,
                            MolForceTable &forces,
                            double scale_force=1) const;
                            
    void calculateAngleForce(const InternalPotential::Molecule &molecule,
                             MolForceTable &forces,
                             double scale_force=1) const;
                             
    void calculateDihedralForce(const InternalPotential::Molecule &molecule,
                                MolForceTable &forces,
                                double scale_force=1) const;

    void calculateImproperForce(const InternalPotential::Molecule &molecule,
                                MolForceTable &forces,
                                double scale_force=1) const;
                                
    void calculateUBForce(const InternalPotential::Molecule &molecule,
                          MolForceTable &forces,
                          double scale_force=1) const;

    void calculateSSForce(const InternalPotential::Molecule &molecule,
                          MolForceTable &forces,
                          double scale_force=1) const;

    void calculateSBForce(const InternalPotential::Molecule &molecule,
                          MolForceTable &forces,
                          double scale_force=1) const;

    void calculateBBForce(const InternalPotential::Molecule &molecule,
                          MolForceTable &forces,
                          double scale_force=1) const;

    void calculateSBTForce(const InternalPotential::Molecule &molecule,
                           MolForceTable &forces,
                           double scale_force=1) const;
                          
    void calculateForce(const InternalPotential::Molecule &molecule,
                        MolForceTable &forces,
                        const Symbol &symbol,
                        const Components &components,
                        double scale_force=1) const;

    bool isstrict;

private:
    void calculatePhysicalEnergy(const GroupInternalParameters &group_params,
                                 const CoordGroup *cgroup_array,
                                 InternalPotential::Energy &energy,
                                 double scale_energy) const;

    void calculateNonPhysicalEnergy(const GroupInternalParameters &group_params,
                                    const CoordGroup *cgroup_array,
                                    InternalPotential::Energy &energy,
                                    double scale_energy) const;

    void calculateCrossEnergy(const GroupInternalParameters &group_params,
                              const CoordGroup *cgroup_array,
                              InternalPotential::Energy &energy,
                              double scale_energy) const;

    static ParameterNames param_names;
    static InternalSymbols internal_symbols;
};

/** This is a forcefield that calculates the energies and forces
    caused by molecular mechanics style internal intramolecular
    potentials, e.g. bond, angle, dihedral, urey bradley terms
    
    @author Christopher Woods
*/  
class SIREMM_EXPORT InternalFF 
            : public SireBase::ConcreteProperty<InternalFF,G1FF>,
              public FF3D,
              protected InternalPotential
{

friend QDataStream& ::operator<<(QDataStream&, const InternalFF&);
friend QDataStream& ::operator>>(QDataStream&, InternalFF&);

public:
    typedef InternalPotential::Components Components;
    typedef InternalPotential::ParameterNames ParameterNames;

    InternalFF();
    InternalFF(const QString &name);
    
    InternalFF(const InternalFF &other);
    
    ~InternalFF();
    
    static const char* typeName();
    
    InternalFF* clone() const;
    
    InternalFF& operator=(const InternalFF &other);
    
    bool operator==(const InternalFF &other) const;
    bool operator!=(const InternalFF &other) const;

    const InternalSymbols& symbols() const;
    const ParameterNames& parameters() const;
    const Components& components() const;

    bool setStrict(bool isstrict);
    bool isStrict() const;

    bool setProperty(const QString &name, const Property &property);
    const Property& property(const QString &name) const;
    bool containsProperty(const QString &name) const;
    const Properties& properties() const;
    
    void force(ForceTable &forcetable, double scale_force=1);
    
    void force(ForceTable &forcetable, const Symbol &symbol,
               double scale_force=1);

    void mustNowRecalculateFromScratch();    
               
protected:
    const FFComponent& _pvt_components() const;

    void recalculateEnergy();

    void _pvt_added(const PartialMolecule &mol, const PropertyMap &map);

    void _pvt_removed(const PartialMolecule &mol);

    void _pvt_changed(const SireMol::Molecule &molecule);
    void _pvt_changed(const QList<SireMol::Molecule> &molecules);
    
    void _pvt_removedAll();
        
    bool _pvt_wouldChangeProperties(MolNum molnum, 
                                    const PropertyMap &map) const;

    void _pvt_updateName();

private:
    typedef InternalPotential::Molecule Molecule;
    typedef InternalPotential::Molecules Molecules;
    typedef InternalPotential::ChangedMolecule ChangedMolecule;
    
    bool recordingChanges() const;
    void recordChange(const ChangedMolecule &change);
    
    /** All of the molecules currently in this forcefield */
    Molecules mols;

    /** The list of molecules that have changed since the last evaluation.
        While ffmols only contains the newest version of the molecule,
        this list contains both the newest version, and the version of the
        molecule at the last energy evaluation. */
    QHash<MolNum,ChangedMolecule> changed_mols;

    /** The energy components available for this forcefield */
    Components ffcomponents;
    
    /** The properties of this forcefield */
    Properties props;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

////////
//////// Inline functions of InternalFF
////////

/** Return all of the symbols used in the internal energy functions */
inline const InternalSymbols& InternalFF::symbols() const
{
    return InternalPotential::symbols();
}

/** Return the names of all of the properties used to store the 
    parameters for this potential */
inline const InternalFF::ParameterNames& InternalFF::parameters() const
{
    return InternalPotential::parameters();
}

/** Return all of the symbols representing the components
    of this forcefield */
inline const InternalFF::Components& InternalFF::components() const
{
    return ffcomponents;
}

/** Return whether or not this strictly include terms that
    involve *only* selected atoms. Otherwise this includes
    terms that involve at least one selected atom */
inline bool InternalFF::isStrict() const
{
    return isstrict;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireMM::InternalFF );

SIRE_EXPOSE_CLASS( SireMM::BondParameterName )
SIRE_EXPOSE_CLASS( SireMM::AngleParameterName )
SIRE_EXPOSE_CLASS( SireMM::DihedralParameterName )
SIRE_EXPOSE_CLASS( SireMM::ImproperParameterName )
SIRE_EXPOSE_CLASS( SireMM::UreyBradleyParameterName )
SIRE_EXPOSE_CLASS( SireMM::StretchStretchParameterName )
SIRE_EXPOSE_CLASS( SireMM::StretchBendParameterName )
SIRE_EXPOSE_CLASS( SireMM::BendBendParameterName )
SIRE_EXPOSE_CLASS( SireMM::StretchBendTorsionParameterName )
SIRE_EXPOSE_CLASS( SireMM::InternalParameterNames )
SIRE_EXPOSE_CLASS( SireMM::InternalParameterNames3D )
SIRE_EXPOSE_CLASS( SireMM::InternalFF )

#endif