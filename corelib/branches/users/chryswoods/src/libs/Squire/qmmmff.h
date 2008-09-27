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

#ifndef SQUIRE_QMMMFF_H
#define SQUIRE_QMMMFF_H

#include "qmmmpotential.h"
#include "qmpotential.h"

#include "SireFF/g2ff.h"
#include "SireFF/ff3d.h"

#include "SireMM/coulombpotential.h"

SIRE_BEGIN_HEADER

namespace Squire
{
class QMMMFF;
}

QDataStream& operator<<(QDataStream&, const Squire::QMMMFF&);
QDataStream& operator>>(QDataStream&, Squire::QMMMFF&);

namespace Squire
{

using SireBase::Property;
using SireBase::Properties;

using SireFF::ForceTable;

using SireCAS::Symbol;

/** This class provides the default names of the element, 
    charge and coordinates properties required by the
    QMMMPotential
*/
class SQUIRE_EXPORT ChargeElementParameterNames3D
           : public ElementParameterName,
             public SireMM::ChargeParameterName,
             public SireFF::detail::Coords3DParameterName
{
public:
    ChargeElementParameterNames3D() 
          : ElementParameterName(), SireMM::ChargeParameterName(),
            SireFF::detail::Coords3DParameterName()
    {}
    
    ~ChargeElementParameterNames3D()
    {}
};

/** This is a forcefield that implements QM/MM using electrostatic
    embedding. This forcefield calculates the energy of a QM molecule
    that is surrounded by a collection of MM point charges
    
    @author Christopher Woods
*/
class SQUIRE_EXPORT QMMMFF
          : public SireBase::ConcreteProperty<QMMMFF,SireFF::G2FF>,
            public SireFF::FF3D,
            protected QMMMPotential<QMPotential,SireMM::InterCoulombPotential>
{

friend QDataStream& ::operator<<(QDataStream&, const QMMMFF&);
friend QDataStream& ::operator>>(QDataStream&, QMMMFF&);

public:
    typedef ChargeElementParameterNames3D Parameters;
    typedef QMPotential::Components Components;

    QMMMFF();
    QMMMFF(const QString &name);
    
    QMMMFF(const QMMMFF &other);
    
    ~QMMMFF();
    
    static const char* typeName()
    {
        return QMetaType::typeName( qMetaTypeId<QMMMFF>() );
    }
    
    const char* what() const
    {
        return QMMMFF::typeName();
    }
    
    QMMMFF& operator=(const QMMMFF &other);
    
    bool operator==(const QMMMFF &other) const;
    bool operator!=(const QMMMFF &other) const;
    
    const Components& components() const;

    Parameters parameters() const
    {
        return Parameters();
    }

    bool setProperty(const QString &name, const Property &property);
    const Property& property(const QString &name) const;
    bool containsProperty(const QString &name) const;
    const Properties& properties() const;

    void mustNowRecalculateFromScratch();    

    void force(ForceTable &forcetable, double scale_force=1);
    
    void force(ForceTable &forcetable, const Symbol &symbol,
               double scale_force=1);

protected:

    ////
    //// Virtual functions from SireMM::InterCoulombPotential and QMPotential
    ////
    
    void changedPotential();

    ////
    //// Virtual functions from SireFF::FF
    ////
    
    const Components& _pvt_components() const;
    
    void recalculateEnergy();
    
    void _pvt_restore(const SireFF::ForceField &ffield);
    
    ////
    //// Virtual functions from SireFF::G2FF
    ////

    void _pvt_added(quint32 group_id, 
                    const SireMol::PartialMolecule &mol, 
                    const SireBase::PropertyMap &map);
                    
    void _pvt_removed(quint32 group_id, 
                      const SireMol::PartialMolecule &mol);
                      
    void _pvt_changed(quint32 group_id, const SireMol::Molecule &mol);
    
    void _pvt_changed(quint32 group_id, const QList<SireMol::Molecule> &mols);
    
    void _pvt_removedAll(quint32 group_id);
    
    bool _pvt_wouldChangeProperties(quint32 group_id, 
                                    SireMol::MolNum molnum, 
                                    const SireBase::PropertyMap &map) const;
};

} // end of namespace Squire

Q_DECLARE_METATYPE( Squire::QMMMFF )

SIRE_EXPOSE_CLASS( Squire::QMMMFF )

SIRE_END_HEADER

#endif
