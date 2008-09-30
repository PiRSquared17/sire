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

#include "qmmmelecembedpotential.h"
#include "qmprogram.h"

#include "latticecharges.h"

#include "SireBase/errors.h"

#include "SireStream/datastream.h"

using namespace Squire;
using namespace SireMM;
using namespace SireStream;

template class QMMMPotential<QMPotential,InterCoulombPotential>;

static const RegisterMetaType<QMMMElecEmbedPotential> r_qmmm( MAGIC_ONLY,
                                                "Squire::QMMMElecEmbedPotential" );
                                                
/** Serialise to a binary datastream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds,
                                      const QMMMElecEmbedPotential &qmmm)
{
    writeHeader(ds, r_qmmm, 1);
    
    ds << static_cast<const QMMMPotential<QMPotential,InterCoulombPotential>&>(qmmm);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds,
                                      QMMMElecEmbedPotential &qmmm)
{
    VersionID v = readHeader(ds, r_qmmm);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
    
        ds >> static_cast<QMMMPotential<QMPotential,InterCoulombPotential>&>(qmmm);

        qmmm.mergeProperties();
    }
    else
        throw version_error(v, "1", r_qmmm, CODELOC);
    
    return ds;
}

/** Merge the properties of the QM and MM forcefields */
void QMMMElecEmbedPotential::mergeProperties()
{
    //we are only interested in the space, qm program and switching function
    //properties
    props = Properties();
    
    props.setProperty("space", this->space());
    props.setProperty("switchingFunction", this->switchingFunction());
    props.setProperty("quantum program", this->quantumProgram());
}

/** Constructor */
QMMMElecEmbedPotential::QMMMElecEmbedPotential()
                       : QMMMPotential<QMPotential,InterCoulombPotential>()
{
    this->mergeProperties();
}

/** Copy constructor */
QMMMElecEmbedPotential::QMMMElecEmbedPotential(const QMMMElecEmbedPotential &other)
                       : QMMMPotential<QMPotential,InterCoulombPotential>(other),
                         props(other.props)
{}

/** Destructor */
QMMMElecEmbedPotential::~QMMMElecEmbedPotential()
{}

/** Copy assignment operator */
QMMMElecEmbedPotential& 
QMMMElecEmbedPotential::operator=(const QMMMElecEmbedPotential &other)
{
    QMMMPotential<QMPotential,InterCoulombPotential>::operator=(other);
    props = other.props;
    
    return *this;
}

/** Return the space within which the molecules in this potential exist */
const SpaceBase& QMMMElecEmbedPotential::space() const
{
    return QMPotential::space();
}

/** Return the switching function that is used to implement the non-bonded
    cutoff */
const SwitchFunc& QMMMElecEmbedPotential::switchingFunction() const
{
    return MMPotential::switchingFunction();
}

/** Return the handle to the quantum chemical program that is used 
    by this potential to calculate the QM energies and forces */
const QMProg& QMMMElecEmbedPotential::quantumProgram() const
{
    return QMPotential::quantumProgram();
}

/** Set the space within which all of the molecules in this potential
    will exist. This returns whether or not this changes the
    potential. */
bool QMMMElecEmbedPotential::setSpace(const SpaceBase &space)
{
    if (QMPotential::setSpace(space))
    {
        this->mergeProperties();
        return true;
    }
    else
        return false;
}

/** Set the switching function that will be used to implement the 
    non-bonded cutoff in the QM/MM interface */
bool QMMMElecEmbedPotential::setSwitchingFunction(const SwitchFunc &switchfunc)
{
    if (MMPotential::setSwitchingFunction(switchfunc))
    {
        this->mergeProperties();
        return true;
    }
    else
        return false;
}

/** Set the handle to the quantum chemical program that will be
    used by this potential to calculate the QM energies and forces.
    This returns whether or not this changes this potential */
bool QMMMElecEmbedPotential::setQuantumProgram(const QMProg &program)
{
    if (QMPotential::setQuantumProgram(program))
    {
        this->mergeProperties();
        return true;
    }
    else
        return false;
}

/** Set the property called 'name' to the value 'value'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
bool QMMMElecEmbedPotential::setProperty(const QString &name, const Property &value)
{
    if (not this->containsProperty(name))
        throw SireBase::missing_property( QObject::tr(
            "There is no property called \"%1\" in the potential %2. "
            "Available properties are [ %3 ].")
                .arg(name).arg(this->what())
                .arg(Sire::toString(props.propertyKeys())),
                    CODELOC );

    if (QMPotential::containsProperty(name))
    {
        if (QMPotential::setProperty(name, value))
        {
            this->mergeProperties();
            return true;
        }
        else
            return false;
    }
    else
    {
        if (MMPotential::setProperty(name, value))
        {
            this->mergeProperties();
            return true;
        }
        else
            return false;
    }
}

/** Return the property called 'name'

    \throw SireBase::missing_property
*/
const Property& QMMMElecEmbedPotential::property(const QString &name) const
{
    return props.property(name);
}

/** Return whether or not this potential contains the property at name 'name' */
bool QMMMElecEmbedPotential::containsProperty(const QString &name) const
{
    return props.hasProperty(name);
}

/** Return all of the properties of this potential */
const Properties& QMMMElecEmbedPotential::properties() const
{
    return props;
}

/** This converts the MM molecules in 'mmmols' into a set of lattice charges
    that surround the QM molecules in 'qmmols' */
LatticeCharges QMMMElecEmbedPotential::getLatticeCharges(const QMMolecules &qmmols,
                                                         const MMMolecules &mmmols) const
{
    return LatticeCharges();
}

/** Calculate the QM forces on the molecules in 'molecules' and add them 
    onto the forces in the passed force table (scaled by the optional
    scaling constant 'scale_force') */
void QMMMElecEmbedPotential::calculateForce(const QMMolecules &qmmols, 
                                            const MMMolecules &mmmols,
                                            ForceTable &forcetable, 
                                            double scale_force) const
{
    throw SireError::incomplete_code( QObject::tr(
        "Need to write the code to calculate QM forces!!!"), CODELOC );
}
                    
/** Calculate the QM forces on the molecules in 'molecules' and add them 
    onto the forces in the passed force table (scaled by the optional
    scaling constant 'scale_force') */
void QMMMElecEmbedPotential::calculateForce(const QMMolecules &qmmols,
                                            const MMMolecules &mmmols, 
                                            ForceTable &forcetable,
                                            const Symbol &symbol, 
                                            const Components &components,
                                            double scale_force) const
{
    throw SireError::incomplete_code( QObject::tr(
        "Need to write the code to calculate QM forces!!!"), CODELOC );
}

/** Calculate the QM energy of the molecules in 'qmmols' in the electrostatic
    field of the molecules in 'mmmols' and add it on to the energy 'nrg', optionally
    scaled by the scaling constant 'scale_energy' */
void QMMMElecEmbedPotential::calculateEnergy(const QMMolecules &qmmols, 
                                             const MMMolecules &mmmols,
                                             Energy &nrg, double scale_energy) const
{
    if (scale_energy == 0)
        return;

    //map all of the molecules so that they are in this space
    QMMolecules mapped_qmmols = QMPotential::mapIntoSpace(qmmols);

    LatticeCharges charges = this->getLatticeCharges(mapped_qmmols, mmmols);

    double qmnrg = this->quantumProgram().calculateEnergy(mapped_qmmols, charges);
    
    nrg += Energy(scale_energy * qmnrg);
}

/** Return the contents of the QM program command file that will be used
    to calculate the QM energy of the molecules in 'qmmols' surrounded
    in the field of the molecules in 'mmols' */
QString QMMMElecEmbedPotential::energyCommandFile(
                          const QMMMElecEmbedPotential::QMMolecules &qmmols,
                          const QMMMElecEmbedPotential::MMMolecules &mmmols) const
{
    //map all of the molecules so that they are in this space
    QMMolecules mapped_qmmols = QMPotential::mapIntoSpace(qmmols);

    LatticeCharges charges = this->getLatticeCharges(mapped_qmmols, mmmols);

    return this->quantumProgram().energyCommandFile(mapped_qmmols, charges);
}

/** Return the contents of the QM program command file that will be used
    to calculate the QM forces on the molecules 'qmmols', and on the 
    molecules 'mmols', where 'qmmols' are the QM molecules that are
    in the field of the MM molecules 'mmmols' */
QString QMMMElecEmbedPotential::forceCommandFile(
                          const QMMMElecEmbedPotential::QMMolecules &qmmols,
                          const QMMMElecEmbedPotential::MMMolecules &mmmols) const
{
    //map all of the molecules so that they are in this space
    QMMolecules mapped_qmmols = QMPotential::mapIntoSpace(qmmols);

    LatticeCharges charges = this->getLatticeCharges(mapped_qmmols, mmmols);

    return this->quantumProgram().forceCommandFile(mapped_qmmols, charges);
}