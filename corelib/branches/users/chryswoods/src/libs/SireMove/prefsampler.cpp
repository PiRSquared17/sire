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

#include "prefsampler.h"

#include "SireMol/atomselection.h"
#include "SireMol/molidentifier.h"

#include "SireBase/majorminorversion.h"

#include "SireSystem/system.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMove;
using namespace SireMol;
using namespace SireBase;
using namespace SireMaths;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<PrefSampler> r_prefsampler;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const PrefSampler &prefsampler)
{
    writeHeader(ds, r_prefsampler, 1);
    
    SharedDataStream sds(ds);
    
    sds << prefsampler.focal_molecule << prefsampler.focal_point 
        << prefsampler.coords_property
        << prefsampler.space_property
        << prefsampler.sampling_constant
        << prefsampler.current_space
        << static_cast<const Sampler&>(prefsampler);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, PrefSampler &prefsampler)
{
    VersionID v = readHeader(ds, r_prefsampler);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> prefsampler.focal_molecule >> prefsampler.focal_point
            >> prefsampler.coords_property
            >> prefsampler.space_property
            >> prefsampler.sampling_constant
            >> prefsampler.current_space
            >> static_cast<Sampler&>(prefsampler);
            
        prefsampler.is_dirty = true;
    }
    else
        throw version_error( v, "1", r_prefsampler, CODELOC );
        
    return ds;
}

Q_GLOBAL_STATIC( QMutex, getMutex );

/** Completely recalculate the weights from scratch */
void PrefSampler::recalculateWeights()
{
    if (not is_dirty)
        return;
    
    //we need to protect access to this function as we break copy-on-write
    //by calling this function from a const-member function
    QMutexLocker lkr( getMutex() );
    
    if (not is_dirty)
        //someone beat us to the recalculation!
        return;
    
    PropertyMap map;
    map.set("coordinates", coords_property);

    //calculate the new center
    if (this->usingFocalMolecule())
    {
        focal_point = focal_molecule.evaluate().centerOfGeometry(map);
    }
        
    //recalculate the weights...
    const MoleculeGroup &molgroup = this->group();
    const QVector< tuple<MolNum,Index> > &viewindicies = molgroup.molViewIndicies();
    int nviews = viewindicies.count();
    
    molweights = QVector<double>( nviews );
    molweights.squeeze();
    
    const tuple<MolNum,Index> *viewindicies_array = viewindicies.constData();
    double *molweights_array = molweights.data();
    
    sum_of_weights = 0;
    max_weight = 0;
    
    for (int i=0; i<nviews; ++i)
    {
        const tuple<MolNum,Index> &viewidx = viewindicies_array[i];
        
        const ViewsOfMol &mol = molgroup[viewidx.get<0>()];
            
        //calculate the distance from the focal point
        Vector new_center = mol.at(viewidx.get<1>()).evaluate()
                               .centerOfGeometry(map);
                               
        double dist2 = current_space->calcDist2(new_center, focal_point);
        
        //weight = 1 / (dist^2 + k)
        double invweight = dist2 + sampling_constant;
        
        if (invweight == 0)
            //default to '1' for zero invweight
            molweights_array[i] = 1;
        else
            molweights_array[i] = 1.0 / invweight;
            
        sum_of_weights += molweights_array[i];
        max_weight = qMax(max_weight, molweights_array[i]);
    }
    
    is_dirty = false;
}

/** Construct a sampler that prefers the molecules near the origin */
PrefSampler::PrefSampler()
            : ConcreteProperty<PrefSampler,Sampler>(),
              coords_property("coordinates"),
              space_property("space"),
              sampling_constant(0), 
              sum_of_weights(0), max_weight(0),
              is_dirty(true)
{}

/** Construct a sampler that prefers the molecules near the origin 
    using the supplied preferential sampling constant 'k' */
PrefSampler::PrefSampler(Area k)
            : ConcreteProperty<PrefSampler,Sampler>(),
              coords_property("coordinates"),
              space_property("space"),
              sampling_constant(k), 
              sum_of_weights(0), max_weight(0),
              is_dirty(true)
{}

/** Construct a sampler that prefers the molecules that are closest
    to the point 'point' */
PrefSampler::PrefSampler(const Vector &point)
            : ConcreteProperty<PrefSampler,Sampler>(),
              focal_point(point),
              coords_property("coordinates"),
              space_property("space"),
              sampling_constant(0),
              sum_of_weights(0), max_weight(0),
              is_dirty(true)
{}

/** Construct a sampler that prefers the molecules that are closest
    to the point 'point' using the supplied preferential
    sampling constant 'k' */
PrefSampler::PrefSampler(const Vector &point, Area k)
            : ConcreteProperty<PrefSampler,Sampler>(),
              focal_point(point),
              coords_property("coordinates"),
              space_property("space"),
              sampling_constant(k),
              sum_of_weights(0), max_weight(0),
              is_dirty(true)
{}

/** Construct a sampler that samples molecule views from the group 'molgroup',
    with a preference to sample molecules that are closest to the 
    point 'point' */
PrefSampler::PrefSampler(const Vector &point, const MoleculeGroup &molgroup)
            : ConcreteProperty<PrefSampler,Sampler>(molgroup),
              focal_point(point), 
              coords_property("coordinates"),
              space_property("space"),
              sampling_constant(0),
              sum_of_weights(0), max_weight(0),
              is_dirty(true)
{}
            
/** Construct a sampler that samples molecule views from the group 'molgroup',
    with a preference to sample molecules that are closest to the 
    point 'point' using the supplied preferential
    sampling constant 'k' */
PrefSampler::PrefSampler(const Vector &point, const MoleculeGroup &molgroup,
                         Area k)
            : ConcreteProperty<PrefSampler,Sampler>(molgroup),
              focal_point(point), 
              coords_property("coordinates"),
              space_property("space"),
              sampling_constant(k),
              sum_of_weights(0), max_weight(0),
              is_dirty(true)
{}

/** Construct a sampler that prefers the molecules that are closest to 
    the view 'molview' */
PrefSampler::PrefSampler(const MoleculeView &molview)
            : ConcreteProperty<PrefSampler,Sampler>(),
              focal_molecule(molview), 
              coords_property("coordinates"),
              sampling_constant(0),
              is_dirty(true)
{}

/** Construct a sampler that prefers the molecules that are closest to 
    the view 'molview' using the supplied preferential
    sampling constant 'k' */
PrefSampler::PrefSampler(const MoleculeView &molview, Area k)
            : ConcreteProperty<PrefSampler,Sampler>(),
              focal_molecule(molview), 
              coords_property("coordinates"),
              sampling_constant(k),
              is_dirty(true)
{}

/** Construct a sampler that samples the molecules in 'molgroup', and
    that prefers molecules that are closest to the view 'molview' */
PrefSampler::PrefSampler(const MoleculeView &molview, const MoleculeGroup &molgroup)
            : ConcreteProperty<PrefSampler,Sampler>(molgroup),
              focal_molecule(molview),
              coords_property("coordinates"),
              space_property("space"),
              sampling_constant(0),
              sum_of_weights(0), max_weight(0),
              is_dirty(true)
{}

/** Construct a sampler that samples the molecules in 'molgroup', and
    that prefers molecules that are closest to the view 'molview' 
    using the optionally supplied preferential sampling constant 'k' */
PrefSampler::PrefSampler(const MoleculeView &molview, const MoleculeGroup &molgroup,
                         Area k)
            : ConcreteProperty<PrefSampler,Sampler>(molgroup),
              focal_molecule(molview),
              coords_property("coordinates"),
              space_property("space"),
              sampling_constant(k),
              sum_of_weights(0), max_weight(0),
              is_dirty(true)
{}

/** Copy constructor */
PrefSampler::PrefSampler(const PrefSampler &other)
            : ConcreteProperty<PrefSampler,Sampler>(other),
              focal_molecule(other.focal_molecule),
              focal_point(other.focal_point),
              coords_property(other.coords_property),
              space_property("space"),
              sampling_constant(other.sampling_constant),
              sum_of_weights(other.sum_of_weights), 
              max_weight(other.max_weight),
              molweights(other.molweights),
              current_space(other.current_space),
              is_dirty(other.is_dirty)
{}

/** Destructor */
PrefSampler::~PrefSampler()
{}

/** Copy assignment operator */
PrefSampler& PrefSampler::operator=(const PrefSampler &other)
{
    if (this != &other)
    {
        focal_molecule = other.focal_molecule;
        focal_point = other.focal_point;
        coords_property = other.coords_property;
        space_property = other.space_property;
        sampling_constant = other.sampling_constant;
        sum_of_weights = other.sum_of_weights;
        max_weight = other.max_weight;
        molweights = other.molweights;
        current_space = other.current_space;
        is_dirty = other.is_dirty;
        
        Sampler::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool PrefSampler::operator==(const PrefSampler &other) const
{
    return focal_molecule == other.focal_molecule and
           focal_point == other.focal_point and
           coords_property == other.coords_property and
           space_property == other.space_property and
           current_space == other.current_space and
           sampling_constant == other.sampling_constant and
           Sampler::operator==(other);
}

/** Comparison operator */
bool PrefSampler::operator!=(const PrefSampler &other) const
{
    return not this->operator==(other);
}

/** Set the molecule group from which molecule view will be sampled */
void PrefSampler::setGroup(const MoleculeGroup &molgroup)
{
    if (molgroup == this->group())
        //nothing to do
        return;

    if (this->usingFocalMolecule())
    {
        //does this group contain the focal molecule?
        MolNum molnum = focal_molecule.data().number();
        
        if (molgroup.contains(molnum))
        {
            if (molgroup[molnum].data().version() != focal_molecule.data().version())
            {
                focal_molecule.update( molgroup[molnum].data() );
                is_dirty = true;
                
                Sampler::setGroup(molgroup);
                return;
            }
        }
    }

    //if this is dirty, then there's nothing we should do now
    //(it all needs recalculating)
    if (is_dirty)
    {
        Sampler::setGroup(molgroup);
        return;
    }

    //has the major version number of the molecule group changed?
    if (molgroup.version().majorVersion() != this->group().version().majorVersion())
    {
        //the actual identity of the views may have changed - it is safest
        //to just recalculate the weights from scratch
        is_dirty = true;
        Sampler::setGroup(molgroup);
        return;
    }

    //ok - only the state of some of the views has changed
    //and the central molecule has not changed
    const QVector< tuple<MolNum,Index> > &viewindicies = molgroup.molViewIndicies();
    int nviews = viewindicies.count();
    
    BOOST_ASSERT( nviews == molweights.count() );
    
    const tuple<MolNum,Index> *viewindicies_array = viewindicies.constData();
    double *molweights_array = molweights.data();
    const MoleculeGroup &current_group = this->group();
    
    PropertyMap map;
    map.set("coordinates", coords_property);
    
    sum_of_weights = 0;
    max_weight = 0;
    
    for (int i=0; i<nviews; ++i)
    {
        const tuple<MolNum,Index> &viewidx = viewindicies_array[i];
        
        const ViewsOfMol &mol = molgroup[viewidx.get<0>()];
        
        if (mol.data().version() == current_group[mol.data().number()].version())
        {
            //the molecule hasn't changed
            sum_of_weights += molweights_array[i];
            max_weight = qMax(max_weight, molweights_array[i]);
            continue;
        }
            
        //the molecule has changed - calculate the new distance from
        //the focal point
        Vector new_center = mol.at(viewidx.get<1>()).evaluate()
                               .centerOfGeometry(map);
                               
        double dist2 = current_space->calcDist2(new_center, focal_point);
        
        //weight = 1 / (dist^2 + k)
        double invweight = dist2 + sampling_constant;
        
        if (invweight == 0)
            //default to '1' for zero invweight
            molweights_array[i] = 1;
        else
            molweights_array[i] = 1.0 / invweight;
            
        sum_of_weights += molweights_array[i];
        max_weight = qMax(max_weight, molweights_array[i]);
    }
        
    Sampler::setGroup(molgroup);
}

/** Update this sampler so that the molecules match the versions
    in the system 'system' */
void PrefSampler::updateFrom(const System &system)
{
    if (this->usingFocalMolecule())
    {
        //has the focal molecule changed?
        MolNum molnum = focal_molecule.data().number();
        
        if (system.contains(molnum))
        {
            if (system[molnum].data().version() != focal_molecule.data().version())
            {
                //yes - the system contains the molecule with a different version
                focal_molecule.update( system[molnum].data() );
                is_dirty = true;
            }
        }
    }

    //get the new space
    const Space &system_space = system.property(space_property).asA<Space>();

    if (current_space != system_space)
    {
        current_space = system_space;
        is_dirty = true;
    }
    
    //now update the MoleculeGroup
    Sampler::updateFrom(system);
}

/** Sample a molecule from the group, and return it, together
    with the probability with which it was chosen */
tuple<PartialMolecule,double> PrefSampler::sample() const
{
    const_cast<PrefSampler*>(this)->recalculateWeights();
    
    //sample the molecule

    //use the von Neumann rejection method to choose a random molecule
    //using the owicki weights
    // (weights = 1 / (dist^2 + k), rejection method - choose random
    //  molecule, then choose random number from 0 to max_prob. If
    //  probability of the molecule <= the random number, then accept
    //  this molecule, else go back to the beginning and try again...

    BOOST_ASSERT(max_weight > 0);

    //pick a random molecule...
    while (true)
    {
        //choose a random molecule ID...
        int i = this->generator().randInt(molweights.count()-1);

        //get the desired probability of choosing this molecule...
        double molweight = molweights.constData()[i];

        //compare the normalised probability against a random number from 0 to max_prob
        if (this->generator().rand(max_weight) <= molweight)
        {
            //test passed - return the molecule and the normalised probability
            //of picking the molecule
            return tuple<PartialMolecule,double>(this->group().viewAt(i),
                                                 molweight / sum_of_weights);
        }
    }

}

/** Return the probability with which the molecule 'molecule' was 
    sampled from this sampler 
    
    \throw SireMol::missing_molecule
    \throw SireError::invalid_index
*/
double PrefSampler::probabilityOf(const PartialMolecule &molecule) const
{
    const_cast<PrefSampler*>(this)->recalculateWeights();
    
    //get the index of this specific group
    int idx = this->group().indexOf(molecule);
    
    if (idx == -1)
    {
        //the molecule was not in the group, so had no chance of being picked
        return 0;
    }
    else
    {
        //return the normalised probability
        return molweights.constData()[idx] / sum_of_weights;
    }
}

/** Set the focal molecule of this sampler */
void PrefSampler::setFocalMolecule(const MoleculeView &molview)
{
    focal_molecule = PartialMolecule(molview);
    is_dirty = true;
}

/** Set the focal point of this sampler */
void PrefSampler::setFocalPoint(const Vector &point)
{
    focal_molecule = PartialMolecule();
    focal_point = point;
    is_dirty = true;
}

/** Set the preferential sampling constant - this is the value
    of 'k' in the biasing algorithm */
void PrefSampler::setSamplingConstant(Area k)
{
    sampling_constant = k;
}

/** Return the preferential sampling constant (k) */
Area PrefSampler::samplingConstant() const
{
    return Area(sampling_constant);
}

/** Return whether we are using a focal molecule (rather
    than a focal point) */
bool PrefSampler::usingFocalMolecule() const
{
    return not focal_molecule.selection().selectedNone();
}

/** Return whether we are using a focal point (rather
    than a focal molecule) */
bool PrefSampler::usingFocalPoint() const
{
    return focal_molecule.selection().selectedNone();
}

/** Return the current focal point (this will be the current center of the 
    molecule view if a focal molecule is being used) */
const Vector& PrefSampler::focalPoint() const
{
    return focal_point;
}

/** Return the current focal molecule (this will be an empty molecule
    if a focal point is being used) */
const PartialMolecule& PrefSampler::focalMolecule() const
{
    return focal_molecule;
}
    
/** Set the property name used to find the coordinates property of
    the molecules */
void PrefSampler::setCoordinatesProperty(const PropertyName &coordsproperty)
{
    coords_property = coordsproperty;
}

/** Set the property used to find the space property of the system */
void PrefSampler::setSpaceProperty(const PropertyName &spaceproperty)
{
    space_property = spaceproperty;
}

/** Return the name of the property used to find the coordinates
    of the molecules */
const PropertyName& PrefSampler::coordinatesProperty() const
{
    return coords_property;
}

/** Return the name of the property used to find the system space */
const PropertyName& PrefSampler::spaceProperty() const
{
    return space_property;
}
