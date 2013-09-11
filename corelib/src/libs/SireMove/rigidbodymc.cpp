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

#include "rigidbodymc.h"
#include "uniformsampler.h"
#include "ensemble.h"

#include "SireSystem/system.h"

#include "SireMol/partialmolecule.h"
#include "SireMol/molecule.h"

#include "SireVol/space.h"

#include "SireMaths/quaternion.h"
#include "SireMaths/vectorproperty.h"

#include "SireUnits/units.h"
#include "SireUnits/temperature.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>
#include <QTime>

using namespace SireMove;
using namespace SireSystem;
using namespace SireMol;
using namespace SireUnits;
using namespace SireVol;
using namespace SireStream;
using namespace SireMaths;

static const RegisterMetaType<RigidBodyMC> r_rbmc;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds,
                                        const RigidBodyMC &rbmc)
{
    writeHeader(ds, r_rbmc, 5);

    SharedDataStream sds(ds);

    sds << rbmc.smplr << rbmc.center_function
        << rbmc.adel << rbmc.rdel
        << rbmc.reflect_center << rbmc.reflect_radius
        << rbmc.reflect_moves
        << rbmc.sync_trans << rbmc.sync_rot << rbmc.common_center
        << static_cast<const MonteCarlo&>(rbmc);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, RigidBodyMC &rbmc)
{
    VersionID v = readHeader(ds, r_rbmc);

    rbmc.center_function = GetCOGPoint();

    if (v == 5)
    {
        SharedDataStream sds(ds);

        sds >> rbmc.smplr >> rbmc.center_function
            >> rbmc.adel >> rbmc.rdel
            >> rbmc.reflect_center >> rbmc.reflect_radius
            >> rbmc.reflect_moves
            >> rbmc.sync_trans >> rbmc.sync_rot
            >> rbmc.common_center
            >> static_cast<MonteCarlo&>(rbmc);
    }
    else if (v == 4)
    {
        SharedDataStream sds(ds);
        
        sds >> rbmc.smplr >> rbmc.center_function
            >> rbmc.adel >> rbmc.rdel
            >> rbmc.sync_trans >> rbmc.sync_rot
            >> rbmc.common_center
            >> static_cast<MonteCarlo&>(rbmc);
    }
    else if (v == 3)
    {
        SharedDataStream sds(ds);
        
        sds >> rbmc.smplr
            >> rbmc.adel >> rbmc.rdel
            >> rbmc.sync_trans >> rbmc.sync_rot
            >> rbmc.common_center
            >> static_cast<MonteCarlo&>(rbmc);
    }
    else if (v == 2)
    {
        SharedDataStream sds(ds);
        
        sds >> rbmc.smplr
            >> rbmc.adel >> rbmc.rdel
            >> rbmc.sync_trans >> rbmc.sync_rot
            >> static_cast<MonteCarlo&>(rbmc);
            
        rbmc.common_center = false;
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> rbmc.smplr
            >> rbmc.adel >> rbmc.rdel
            >> static_cast<MonteCarlo&>(rbmc);
            
        rbmc.sync_trans = false;
        rbmc.sync_rot = false;
        rbmc.common_center = false;
    }
    else
        throw version_error(v, "1-5", r_rbmc, CODELOC);

    return ds;
}

/** Null constructor */
RigidBodyMC::RigidBodyMC(const PropertyMap &map) 
            : ConcreteProperty<RigidBodyMC,MonteCarlo>(map),
              center_function( GetCOGPoint() ),
              adel( 0.15 * angstrom ), rdel( 15 * degrees ),
              reflect_center(0), reflect_radius(0), reflect_moves(false),
              sync_trans(false), sync_rot(false), common_center(false)
{
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
}

/** Construct a move that moves molecules returned by the sampler 'sampler' */
RigidBodyMC::RigidBodyMC(const Sampler &sampler, const PropertyMap &map)
            : ConcreteProperty<RigidBodyMC,MonteCarlo>(map),
              smplr(sampler), center_function( GetCOGPoint() ),
              adel( 0.15 * angstrom ),
              rdel( 15 * degrees ), 
              reflect_center(0), reflect_radius(0), reflect_moves(false),
              sync_trans(false), sync_rot(false), common_center(false)
{
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
    smplr.edit().setGenerator( this->generator() );
}

/** Construct a move that moves molecule views from the molecule group 'molgroup',
    selecting views randomly, with each view having an equal chance of
    being chosen */
RigidBodyMC::RigidBodyMC(const MoleculeGroup &molgroup, 
                         const PropertyMap &map)
            : ConcreteProperty<RigidBodyMC,MonteCarlo>(map), 
              smplr( UniformSampler(molgroup) ),
              center_function( GetCOGPoint() ),
              adel( 0.15 * angstrom ), rdel( 15 * degrees ),
              reflect_center(0), reflect_radius(0), reflect_moves(false),
              sync_trans(false), sync_rot(false), common_center(false)
{
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
    smplr.edit().setGenerator( this->generator() );
}

/** Copy constructor */
RigidBodyMC::RigidBodyMC(const RigidBodyMC &other)
            : ConcreteProperty<RigidBodyMC,MonteCarlo>(other), 
              smplr(other.smplr), center_function(other.center_function),
              adel(other.adel), rdel(other.rdel),
              reflect_center(other.reflect_center), 
              reflect_radius(other.reflect_radius),
              reflect_moves(other.reflect_moves),
              sync_trans(other.sync_trans), sync_rot(other.sync_rot),
              common_center(other.common_center)
{}

/** Destructor */
RigidBodyMC::~RigidBodyMC()
{}

void RigidBodyMC::_pvt_setTemperature(const Temperature &temperature)
{
    MonteCarlo::setEnsemble( Ensemble::NVT(temperature) );
}

/** Copy assignment operator */
RigidBodyMC& RigidBodyMC::operator=(const RigidBodyMC &other)
{
    if (this != &other)
    {
        smplr = other.smplr;
        center_function = other.center_function;
        adel = other.adel;
        rdel = other.rdel;
        reflect_center = other.reflect_center;
        reflect_radius = other.reflect_radius;
        reflect_moves = other.reflect_moves;
        sync_trans = other.sync_trans;
        sync_rot = other.sync_rot;
        common_center = other.common_center;
        MonteCarlo::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool RigidBodyMC::operator==(const RigidBodyMC &other) const
{
    return smplr == other.smplr and center_function == other.center_function and
           adel == other.adel and rdel == other.rdel and
           reflect_center == other.reflect_center and
           reflect_radius == other.reflect_radius and
           reflect_moves == other.reflect_moves and 
           sync_trans == other.sync_trans and sync_rot == other.sync_rot and
           common_center == other.common_center and
           MonteCarlo::operator==(other);
}

/** Comparison operator */
bool RigidBodyMC::operator!=(const RigidBodyMC &other) const
{
    return not this->operator==(other);
}

/** Return a string representation of this move */
QString RigidBodyMC::toString() const
{
    return QObject::tr("RigidBodyMC( maximumTranslation() = %1 A, "
                       "maximumRotation() = %2 degrees "
                       "nAccepted() = %3 nRejected() = %4 )")
                  .arg(this->maximumTranslation().to(angstrom))
                  .arg(this->maximumRotation().to(degrees))
                  .arg(this->nAccepted())
                  .arg(this->nRejected());
}

/** Set the maximum delta for any translation */
void RigidBodyMC::setMaximumTranslation(Dimension::Length max_translation)
{
    adel = max_translation;
}

/** Set the maximum delta for any rotation */
void RigidBodyMC::setMaximumRotation(Dimension::Angle max_rotation)
{
    rdel = max_rotation;
}

/** Set the function used to get the center of rotation for each molecule */
void RigidBodyMC::setCenterOfRotation(const GetPoint &func)
{
    center_function = func;
}

/** Return the maximum translation for each move */
Dimension::Length RigidBodyMC::maximumTranslation() const
{
    return Dimension::Length(adel);
}

/** Return the maximum rotation for each move */
Dimension::Angle RigidBodyMC::maximumRotation() const
{
    return rdel;
}

/** Return the function used to get the center of rotation of each molecule */
const GetPoint& RigidBodyMC::centerOfRotation() const
{
    return center_function.read();
}

/** Set the sampler (and contained molecule group) that provides
    the random molecules to be moved. This gives the sampler the
    same random number generator that is used by this move */
void RigidBodyMC::setSampler(const Sampler &sampler)
{
    smplr = sampler;
    smplr.edit().setGenerator( this->generator() );
}

/** Set the sampler to be one that selects views at random 
    from the molecule group 'molgroup' (each view has an
    equal chance of being chosen) */
void RigidBodyMC::setSampler(const MoleculeGroup &molgroup)
{
    smplr = UniformSampler(molgroup);
}

/** Return the sampler that is used to draw random molecules */
const Sampler& RigidBodyMC::sampler() const
{
    return smplr;
}

/** Return the molecule group from which random molecule views 
    are drawn */
const MoleculeGroup& RigidBodyMC::moleculeGroup() const
{
    return smplr->group();
}

/** Set the random number generator used by this move */
void RigidBodyMC::setGenerator(const RanGenerator &rangenerator)
{
    MonteCarlo::setGenerator(rangenerator);
    smplr.edit().setGenerator(this->generator());
}

/** Set whether or not to synchronise translation of all of the views */
void RigidBodyMC::setSynchronisedTranslation(bool on)
{
    sync_trans = on;
}

/** Set whether or not to synchronise rotation of all of the views */
void RigidBodyMC::setSynchronisedRotation(bool on)
{
    sync_rot = on;
}

/** Set whether or not to use the same rotation center for all
    synchronised molecules */
void RigidBodyMC::setSharedRotationCenter(bool on)
{
    common_center = on;
}

/** Turn on rigid body move reflections. This makes sure that only
    molecules within the specified sphere can be moved, and any moves
    are reflected so that these molecules will always remain within 
    the sphere */
void RigidBodyMC::setReflectionSphere(Vector sphere_center,
                                      SireUnits::Dimension::Length sphere_radius)
{
    reflect_moves = true;
    reflect_center = sphere_center;
    reflect_radius = sphere_radius.value();

    if (reflect_radius < 0.01)
    {
        reflect_moves = false;
        reflect_center = Vector(0);
        reflect_radius = 0;
    }
}

/** Return whether or not these moves use a reflection sphere */
bool RigidBodyMC::usesReflectionMoves() const
{
    return reflect_moves;
}

/** Return the center of the reflection sphere. Returns a null vector
    if a reflection sphere is not being used */
Vector RigidBodyMC::reflectionSphereCenter() const
{
    return reflect_center;
}

/** Return the radius of the reflection sphere. This returns zero
    if the reflection sphere is not being used */
SireUnits::Dimension::Length RigidBodyMC::reflectionSphereRadius() const
{
    return SireUnits::Dimension::Length( reflect_radius );
}

/** Return whether or not translation of all molecules is synchronised */
bool RigidBodyMC::synchronisedTranslation() const
{
    return sync_trans;
}

/** Return whether or not rotation of all molecules is synchronised */
bool RigidBodyMC::synchronisedRotation() const
{
    return sync_rot;
}

/** Return whether or not synchronised rotation uses the same
    center of rotation for all molecules */
bool RigidBodyMC::sharedRotationCenter() const
{
    return common_center;
}

/** This internal function is used to actually move the molecule(s) */
void RigidBodyMC::performMove(System &system,
                              double &old_bias, double &new_bias,
                              const PropertyMap &map)
{
    const PropertyName &center_property = map["center"];

    //update the sampler with the latest version of the molecules
    smplr.edit().updateFrom(system);

    //get the random amounts by which to translate and
    //rotate the molecule(s)
    Vector delta = generator().vectorOnSphere(adel);

    const Quaternion rotdelta( rdel * generator().rand(),
                               generator().vectorOnSphere() );

    if ( (not sync_trans) and (not sync_rot) )
    {
        //randomly select a molecule to move
        tuple<PartialMolecule,double> mol_and_bias = smplr.read().sample();

        const PartialMolecule &oldmol = mol_and_bias.get<0>();
        old_bias = mol_and_bias.get<1>();
        
        if (oldmol.isEmpty())
        {
            qDebug() << "Sampler returned an empty molecule in RigidBodyMC" << this->toString()
                     << this->moleculeGroup().toString()
                     << this->moleculeGroup().nMolecules() << smplr.read().toString();
            return;
        }
        
        const bool has_center_property = (oldmol.selectedAll() and
                                          oldmol.hasProperty(center_property));

        //move the molecule
        PartialMolecule newmol;
        
        if (has_center_property)
        {
            newmol = oldmol.move()
                           .rotate(rotdelta,
                                   oldmol.property(center_property).asA<VectorProperty>(),
                                   map)
                           .translate(delta, map)
                           .commit();
        }
        else
        {
            newmol = oldmol.move()
                           .rotate(rotdelta,
                                   center_function.read()(oldmol,map),
                                   map)
                           .translate(delta, map)
                           .commit();
        }

        //if we are reflecting moves in a sphere, then check that this move
        //won't take us out of the sphere.
        if (reflect_moves)
        {
            //moves are constrained into a sphere of radius reflect_radius
            //around reflect_center. If the center of geometry moves outside
            //the sphere, then the molecule will bounce off the edge of 
            //the sphere and back into the sphere volume

            Vector old_center;
            
            if (has_center_property)
            {
                old_center = oldmol.property(center_property).asA<VectorProperty>();
            }
            else
            {
                old_center = oldmol.evaluate().center();
            }

            if ( (old_center-reflect_center).length() > reflect_radius )
            {
                //the molecule is already outside the sphere, so cannot be moved
                old_bias = 1;
                new_bias = 1;
                qDebug() << "HOW IS THE MOLECULE OUTSIDE THE SPHERE?";
                qDebug() << (old_center-reflect_center).length() << reflect_radius;
                return;
            }

            Vector new_center;
            
            if (has_center_property)
            {
                new_center = newmol.property(center_property).asA<VectorProperty>();
            }
            else
            {
                new_center = newmol.evaluate().center();
            }

            double dist = (new_center - reflect_center).length();

            if (dist > reflect_radius)
            {
                delta = new_center - old_center;
            
                //this would move the molecule out of the sphere. We need
                //to work out where the molecule would intersect the surface
                //of the sphere, and then reflect the molecule back inside

                //first, find the intersection of the delta vector with the sphere.
                // The delta vector has origin at O, direction D. The sphere
                // is at origin C, with radius R
                Vector D = delta.normalise();
                Vector O = old_center;
                Vector C = reflect_center;
                double R2 = reflect_radius*reflect_radius;

                //a point P is on the surface of the sphere if (P-C).(P-C) = R^2
                //this means that the intersection of the vector with the sphere
                //must satisfy ( (O + xD) - C ).( (O + xD) - C ) = R^2
                // This gives;
                // (D.D) x^2 + 2 ( O-C ).D x + (O-C).(O-C) - R^2 = 0
                // which is A x^2 + B x + C = 0, where
                //
                // A = D.D
                // B = 2 (O-C).D
                // C = (O-C).(O-C) - R^2
                //
                // which can be solved using the quadratic formula
                //
                // roots = [-B - sqrt(B^2 - 4AC)] / 2A
                //       = [-B + sqrt(B^2 - 4AC)] / 2A
                //
                // To avoid numerical instability, we use;
                //
                // roots = Q / A and C / Q where
                //
                // Q = [-B + sqrt(B^2 - 4AC)] / 2   if B < 0
                // Q = [-B - sqrt(B^2 - 4AC)] / 2   otherwise

                double QA = Vector::dot(D,D);
                double QB = 2.0 * Vector::dot( O-C, D );
                double QC = Vector::dot(O-C, O-C) - R2;

                double B2_minus_4AC = QB*QB - 4*QA*QC;

                if (B2_minus_4AC < 0)
                {
                    //the move does not intersect with the sphere... weird...
                    old_bias = 1;
                    new_bias = 1;
                    qDebug() << "WEIRD: MOVE DOES NOT INTERSECT WITH SPHERE";
                    return;
                }

                double Q;

                if (QB < 0)
                {
                    Q = (-QB - std::sqrt(B2_minus_4AC)) * 0.5;
                }
                else
                {
                    Q = (-QB + std::sqrt(B2_minus_4AC)) * 0.5;
                }

                double x0 = Q / QA;
                double x1 = QC / Q;

                if (x0 > x1){ qSwap(x0,x1); }

                if (x1 < 0)
                {
                    //the intersection is behind us...
                    qDebug() << "Intersection behind us..." << x1;
                    old_bias = 1.0;
                    new_bias = 1.0;
                    return;
                }

                double x = x0;

                if (x0 < 0){ x = x1; }

                //the intersection point, X, is O + xD
                Vector X = O + x*D;

                //qDebug() << "Reflect at " << X.toString() << (X-C).length() << std::sqrt(R2);

                //ok - now we have the intersection point, the next step is to 
                //get the normal (N) to the sphere at this point, as this defines the
                //reflection plane
                Vector N = (X - C).normalise();

                //We want to reflect the unit vector that intersects with the 
                //sphere about this normal
                Vector X1 = X - D;
                Vector X2 = X1 + 2 * ( (X - Vector::dot(D,N)*N) - X1 );

                //X2 is the reflected vector. Now work out how much we have
                //moved along X1, and then move that same amount along X2
                double dist_x1 = (X-O).length();
                double dist_x2 = delta.length() - dist_x1;

                if (dist_x2 < 0)
                {
                    qDebug() << "WEIRD. NEGATIVE REFLECTION DISTANCE??? " 
                             << dist_x2;
                    
                    old_bias = 1.0;
                    new_bias = 1.0;
                    return;
                }

                //work out where the new center of the molecule should lie
                Vector new_new_center = X + dist_x2*((X2-X).normalise());

                //now translate the molecule to the new position
                newmol = newmol.move().translate(new_new_center-new_center,map).commit();

                new_center = newmol.evaluate().center();
                double dist = (new_center - reflect_center).length();
                
                int check_count = 0;
                
                while (dist > reflect_radius)
                {
                    qDebug() << "MOVED MOLECULE OUTSIDE SPHERE" << dist << reflect_radius;
                    qDebug() 
                        << "FIXING THE PROBLEM (MOSTLY CAUSED BY NUMERICAL IMPRECISION)";
                    
                    //this will be due to a little numerical imprecision
                    newmol = newmol.move().translate( 
                            (1.01*(dist-reflect_radius))
                                * ((reflect_center-new_center).normalise()) ).commit();
                    
                    if (has_center_property)
                    {
                        new_center = newmol.property(center_property).asA<VectorProperty>();
                    }
                    else
                    {
                        new_center = newmol.evaluate().center();
                    }
                    
                    dist = (new_center - reflect_center).length();
                    
                    check_count += 1;
                    
                    if (check_count > 100)
                    {
                        qDebug() << "WARNING: SOMETHING WEIRD GOING ON."
                                 << "SKIPPING THIS MOVE.";
                                 
                        old_bias = 1;
                        new_bias = 1;
                        return;
                    }
                }
            } 
        }

        //update the system with the new coordinates
        system.update(newmol);

        //get the new bias on this molecule
        new_bias = smplr.read().probabilityOf(newmol);
    }
    else if (sync_trans)
    {
        if (sync_rot)
        {
            //translate and rotate all molecules
            old_bias = 1;
            new_bias = 1;

            const Molecules &molecules = smplr.read().group().molecules();

            Molecules new_molecules = molecules;

            if (common_center)
            {
                AABox box;
            
                //rotate all molecules around the same center
                for (Molecules::const_iterator it = molecules.constBegin();
                     it != molecules.constEnd();
                     ++it)
                {
                    if (it->selectedAll() and it->molecule().hasProperty(center_property))
                    {
                        box += it->molecule().property(center_property).asA<VectorProperty>();
                    }
                    else
                    {
                        box += center_function.read()(*it,map);
                    }
                }
                
                for (Molecules::const_iterator it = molecules.constBegin();
                     it != molecules.constEnd();
                     ++it)
                {
                    PartialMolecule newmol = it->move()
                                                .rotate(rotdelta, box.center(), map)
                                                .translate(delta, map)
                                                .commit();
                    
                    new_molecules.update(newmol);
                }            
                
            }
            else
            {
                for (Molecules::const_iterator it = molecules.constBegin();
                     it != molecules.constEnd();
                     ++it)
                {
                    PartialMolecule newmol = it->move()
                                                .rotate(rotdelta,
                                                        center_function.read()(*it,map),
                                                        map)
                                                .translate(delta, map)
                                                .commit();
                    
                    new_molecules.update(newmol);
                }            
            }
            
            system.update(new_molecules);            
        }
        else
        {
            //translate all molecules
            const Molecules &molecules = smplr.read().group().molecules();

            Molecules new_molecules = molecules;

            for (Molecules::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                PartialMolecule newmol = it->move()
                                            .translate(delta, map)
                                            .commit();

                new_molecules.update(newmol);
            }

            system.update(new_molecules);

            //then rotate a single random molecule
            smplr.edit().updateFrom(system);

            tuple<PartialMolecule,double> mol_and_bias = smplr.read().sample();

            const PartialMolecule &oldmol = mol_and_bias.get<0>();
            old_bias = mol_and_bias.get<1>();

            PartialMolecule newmol;
            
            if (oldmol.selectedAll() and oldmol.hasProperty(center_property))
            {
                newmol = oldmol.move()
                               .rotate(rotdelta,
                                       oldmol.property(center_property).asA<VectorProperty>(),
                                       map)
                               .commit();
            }
            else
            {
                newmol = oldmol.move()
                               .rotate(rotdelta,
                                       center_function.read()(oldmol,map),
                                       map)
                               .commit();
            }

            //update the system with the new coordinates
            system.update(newmol);

            //get the new bias on this molecule
            new_bias = smplr.read().probabilityOf(newmol);
        }
    }
    else if (sync_rot)
    {
        //rotate all of the molecules
        const Molecules &molecules = smplr.read().group().molecules();

        Molecules new_molecules = molecules;

        bool perform_translation = true;

        if (common_center)
        {
            //we cannot perform a move with synchronised rotation and
            //individual rotation if a common center is used, as it 
            //would be very hard to work out the probability of the 
            //reverse move... So we will instead have a 50/50 chance
            //of performing rotation only or translation only
            perform_translation = generator().randBool();
            
            if (not perform_translation)
            {
                AABox box;
                
                for (Molecules::const_iterator it = molecules.constBegin();
                     it != molecules.constEnd();
                     ++it)
                {
                    if (it->selectedAll() and it->molecule().hasProperty(center_property))
                    {
                        box += it->molecule().property(center_property).asA<VectorProperty>();
                    }
                    else
                    {
                        box += center_function.read()(*it,map);
                    }
                }
            
                for (Molecules::const_iterator it = molecules.constBegin();
                     it != molecules.constEnd();
                     ++it)
                {
                    PartialMolecule newmol = it->move()
                                                .rotate(rotdelta, box.center(), map)
                                                .commit();

                    new_molecules.update(newmol);
                }
            }
        }
        else
        {
            for (Molecules::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                if (it->selectedAll() and it->molecule().hasProperty(center_property))
                {
                    PartialMolecule newmol = it->move()
                                                .rotate(rotdelta,
                                                        it->molecule().property(center_property)
                                                                .asA<VectorProperty>(),
                                                        map)
                                                .commit();
                    
                    new_molecules.update(newmol);
                }
                else
                {
                    PartialMolecule newmol = it->move()
                                                .rotate(rotdelta,
                                                        center_function.read()(*it,map),
                                                        map)
                                                .commit();
                    
                    new_molecules.update(newmol);
                }
            }
        }

        system.update(new_molecules);

        if (perform_translation)
        {
            //then translate a single random molecule
            smplr.edit().updateFrom(system);

            tuple<PartialMolecule,double> mol_and_bias = smplr.read().sample();

            const PartialMolecule &oldmol = mol_and_bias.get<0>();
            old_bias = mol_and_bias.get<1>();

            PartialMolecule newmol = oldmol.move()
                                           .translate(delta, map)
                                           .commit();

            //update the system with the new coordinates
            system.update(newmol);

            //get the new bias on this molecule
            new_bias = smplr.read().probabilityOf(newmol);
        }
        else
        {
            old_bias = 1;
            new_bias = 1;
        }

    }
}

/** Attempt 'n' rigid body moves of the views of the system 'system' */
void RigidBodyMC::move(System &system, int nmoves, bool record_stats)
{
    if (nmoves <= 0)
        return;

    //save our, and the system's, current state
    RigidBodyMC old_state(*this);

    System old_system_state(system);
    
    try
    {
        const PropertyMap &map = Move::propertyMap();
            
        for (int i=0; i<nmoves; ++i)
        {
            //get the old total energy of the system
            double old_nrg = system.energy( this->energyComponent() );

            //save the old system
            System old_system(system);

            double old_bias = 1;
            double new_bias = 1;

            this->performMove(system, old_bias, new_bias, map);
    
            //calculate the energy of the system
            double new_nrg = system.energy( this->energyComponent() );

            //accept or reject the move based on the change of energy
            //and the biasing factors
            if (not this->test(new_nrg, old_nrg, new_bias, old_bias))
            {
                //the move has been rejected - reset the state
                system = old_system;
            }

            if (record_stats)
            {
                system.collectStats();
            }
        }
    }
    catch(...)
    {
        system = old_system_state;
        this->operator=(old_state);

        throw;
    }
}

const char* RigidBodyMC::typeName()
{
    return QMetaType::typeName( qMetaTypeId<RigidBodyMC>() );
}
