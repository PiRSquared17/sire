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

#include "mover.h"

using namespace SireMol;

/** The default weighting function used by Mover */
static Property default_weightfunc( new RelFromNumber() );

/** Constructor */
MoverBase::MoverBase()
{}

/** Copy constructor */
MoverBase::MoverBase(const MoverBase &other)
          : movable_atoms(other.movable_atoms)
{}

/** Destructor */
MoverBase::~MoverBase()
{}

/** Copy assignment operator */
MoverBase& MoverBase::operator=(const MoverBase &other)
{
    movable_atoms = other.movable_atoms;
    return *this;
}

/** Translate the selected atoms from 'coords' by 'delta'.
    This function assumes that 'selected_atoms' is compatible
    with 'coords' */
void MoverBase::translate(AtomicCoords &coords,
                          const AtomSelection &selected_atoms,
                          const Vector &delta)
{
    if (delta.isZero() or selected_atoms.selectedNone())
        return;

    CoordGroup *coords_array = coords.data();
    int ncg = coords.count();
    
    if (selected_atoms.selectedAll())
    {
        //we are moving everything!
        for (int i=0; i<ncg; ++i)
        {
            coords_array[i] = coords_array[i].edit().translate(delta);
        }
    }
    else if (selected_atoms.selectedAllCutGroups())
    {
        //we are moving all of the CutGroups
        for (CGIdx i(0); i<ncg; ++i)
        {
            if (selected_atoms.selectedAll(i))
            {
                coords_array[i] = coords_array[i].edit().translate(delta);
            }
            else
            {
                QSet<Index> atoms_to_move = selected_atoms.selectedAtoms(i);
                
                coords_array[i] = coords_array[i].edit().translate(atoms_to_move,
                                                                   delta);
            }
        }
    }
    else
    {
        //we are moving only some CutGroups
        QSet<CGIdx> cg_to_move = selected_atoms.selectedCutGroups();
        
        foreach (CGIdx i, cg_to_move)
        {
            if (selected_atoms.selectedAll(i))
            {
                coords_array[i] = coords_array[i].edit().translate(delta);
            }
            else
            {
                QSet<Index> atoms_to_move = selected_atoms.selectedAtoms(i);

                coords_array[i] = coords_array[i].edit().translate(atoms_to_move,
                                                                   delta);
            }
        }
    }
}

/** Rotate the coordinates (in 'coords') of the specified selected
    atoms using the rotation matrix 'rotmat' about the point 'point'.
    This function assumes that coords and selected_atoms are compatible */
void MoverBase::rotate(AtomicCoords &coords,
                       const AtomSelection &selected_atoms,
                       const Matrix &rotmat,
                       const Vector &point)
{
    if (selected_atoms.selectedNone())
        return;
    
    CoordGroup *coords_array = coords.data();
    int ncg = coords.count();
                
    if (selected_atoms.selectedAll())
    {
        //we are rotating everything
        for (int i=0; i<ncg; ++i)
        {
            coords_array[i] = coords_array[i].edit().rotate(rotmat,point);
        }
    }
    else if (selected_atoms.selectedAllCutGroups())
    {
        //we are rotating every CutGroup
        for (CGIdx i(0); i<ncg; ++i)
        {
            if (selected_atoms.selectedAll(i))
            {
                coords_array[i] = coords_array[i].edit().rotate(rotmat,point);
            }
            else
            {
                QSet<Index> atoms_to_move = selected_atoms.selectedAtoms(i);
                
                coords_array[i] = coords_array[i].edit().rotate(atoms_to_move,
                                                                rotmat, delta);
            }
        }
    }
    else
    {
        QSet<CGIdx> cg_to_move = selected_atoms.selectedCutGroups();
        
        foreach (CGIdx i, cg_to_move)
        {
            if (selected_atoms.selectedAll(i))
            {
                coords_array[i] = coords_array[i].edit().rotate(rotmat, point);
            }
            else
            {
                QSet<Index> atoms_to_move = selected_atoms.selectedAtoms(i);
                
                coords_array[i] = coords_array[i].edit().rotate(atoms_to_move,
                                                                rotmat, point);
            }
        }
    }
}

/** Map the selected atoms from 'coords' into the axis set in 'axes'

    This function assumes that coords and selected_atoms are compatible!
*/
void MoverBase::mapInto(AtomicCoords &coords,
                        const AtomSelection &selected_atoms,
                        const AxisSet &axes)
{
    if (selected_atoms.selectedNone())
        return;
        
    CoordGroup *coords_array = coords.data();
    int ncg = coords.count();

    if (selected_atoms.selectedAll())
    {
        //we are moving everything
        for (int i=0; i<ncg; ++i)
        {
            coords_array[i] = coords_array[i].edit().mapInto(axes);
        }
    }
    else if (selected_atoms.selectedAllCutGroups())
    {
        for (CGIdx i(0); i<ncg; ++i)
        {
            if (selected_atoms.selectedAll(i))
            {
                coords_array[i] = coords_array[i].edit().mapInto(axes);
            }
            else
            {
                QSet<Index> atoms_to_move = selected_atoms.selectedAtoms(i);
                
                coords_array[i] = coords_array[i].edit().mapInto(atoms_to_move,
                                                                 axes);
            }
        }
    }
    else
    {
        QSet<CGIdx> cg_to_move = selected_atoms.selectedCutGroups();
        
        foreach (CGIdx i, cg_to_move)
        {
            if (selected_atoms.selectedAll(i))
            {
                coords_array[i] = coords_array[i].edit().mapInto(axes);
            }
            else
            {
                QSet<Index> atoms_to_move = selected_atoms.selectedAtoms(i);
                
                coords_array[i] = coords_array[i].edit().mapInto(atoms_to_move,
                                                                 axes);
            }
        }
    }
}

/** Translate the selected atoms in the molecule whose data is in 'moldata'
    by 'delta', using 'coord_property' to get the coordinates to 
    be translated. This function assumes that selected_atoms
    and moldata are compatible.
    
    \throw SireBase::missing_property
*/
void MoverBase::translate(MoleculeData &moldata 
                          const AtomSelection &selected_atoms,
                          const Vector &delta, 
                          const PropertyMap &map)
{
    //which property contains the coordinates?
    PropertyName coord_property = map["coordinates"];

    //get the current coordinates
    AtomicCoords coords = moldata.property(coord_property);
    
    //translate the coordinates of the selected atoms
    MoverBase::translate(coords, selected_atoms, delta);
    
    //set the new property
    view.setProperty(coord_property, coords);
}

/** Rotate the selected atoms in the molecule whose data
    is in 'moldata' using the rotation matrix 'rotmat' 
    around the point 'point', using 'coord_property'
    to find the property containing the coordinates 
    to be rotated.
    
    This function assumes that moldata and selected_atoms
    are compatible.
    
    \throw SireBase::missing_property
*/
void MoverBase::rotate(MoleculeData &moldata,
                       const AtomSelection &selected_atoms,
                       const Matrix &rotmat,
                       const Vector &point,
                       const PropertyMap &map)
{
    //get the name of the property that contains the coordinates
    PropertyName coord_property = map["coordinates"];

    //get the coordinates to be rotated
    AtomicCoords coords = moldata.property(coord_property);
    
    //rotate the coordinates
    MoverBase::rotate(coords, selected_atoms, rotmat, point);
    
    //set the new property
    moldata.setProperty(coord_property, coords);
}

/** This function maps the selected atoms from their current
    (cartesian) coordinate frame into the coordinate frame
    described by 'axes'. This function assumes that 
    moldata and 'selected_atoms' are compatible
    
    \throw SireBase::missing_property
*/
void MoverBase::mapInto(MoleculeData &moldata,
                        const AtomSelection &selected_atoms,
                        const AxisSet &axes, 
                        const PropertyMap &map)
{
    //get the name of the property that holds the coordinates
    PropertyName coord_property = map["coordinates"];
    
    //get the coordinates to be mapped
    AtomicCoords coords = moldata.property(coord_property);
    
    //map the coordinates
    MoverBase::mapInto(coords, selected_atoms, axes);
    
    //save the new coordinates
    moldata.setProperty(coord_property, coords);
}

/** Translate atoms we are allowed to move from the molecule whose
    data is in 'moldata' by 'delta', finding the coordinates
    using 'coord_property' 
    
    \throw SireBase::missing_property
*/
void MoverBase::translate(MoleculeData &moldata,      
                          const Vector &delta,
                          const PropertyMap &map) const
{
    MoverBase::translate(moldata, movable_atoms, delta, map);
}

/** Rotate the atoms we are allowed to move from the molecule whose
    data is in 'moldata' using the quaternion 'quat' about the 
    point 'point', finding the coordinates using 'coord_property'
    
    \throw SireBase::missing_property
*/
void MoverBase::rotate(MoleculeData &moldata,
                       const Quaternion &quat,
                       const Vector &point,
                       const PropertyMap &map) const
{
    MoverBase::rotate(moldata, movable_atoms, quat.toMatrix(),
                      point, map);
}

/** Rotate the atoms we are allowed to move from the molecule whose
    data is in 'moldata' using the rotation matrix 'rotmat' about the 
    point 'point', finding the coordinates using 'coord_property'
    
    \throw SireBase::missing_property
*/
void MoverBase::rotate(MoleculeData &moldata,
                       const Matrix &rotmat,
                       const Vector &point,
                       const PropertyMap &map) const
{
    MoverBase::rotate(moldata, movable_atoms, rotmat,
                      point, coord_property);
}

/** Change the length of the bond identified by 'bond' by 'delta',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties 
    
    This only moves the movable atoms in this view, and 
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.
    
    The bond is labelled;
    
    atom0--atom1
    
    The molecule is split into two about this bond, i.e.
    atom0 and everything it is bonded to is in group0, while
    atom1 and everything it is bonded to is in group1.
    
    The two groups are then translated along the vector atom0->atom1
    
    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::change(MoleculeData &moldata, const BondID &bond,
                       SireUnits::Dimension::Length delta,
                       const PropertyMap &map) const
{
    if (delta == 0)
        return;

    //get the indicies of the atoms that comprise the bond
    tuple<AtomIdx,AtomIdx> atomidxs = bond.map(moldata.info());

    AtomIdx atom0 = atomidxs.get<0>();
    AtomIdx atom1 = atomidxs.get<1>();

    //get the connectivity property that is used to split
    //the molecule into two parts
    Connectivity connectivity = moldata.property(map["connectivity"]);

    //see if there are any anchors that must be applied to 
    //the section of molecule
    AtomSelection anchors;
    
    if (map.specified("anchors"))
        anchors = moldata.property(map["anchors"]);

    //split the molecule into the two parts that are
    //going to move - the two groups are only able to
    //contain the atoms that are in 'movable_atoms'
    tuple<AtomSelection,AtomSelection> groups = 
                        connectivity.split(atom0, atom1, 
                                           movable_atoms, anchors);

    const AtomSelection &group0 = groups.get<0>();
    const AtomSelection &group1 = groups.get<1>();

    double weight0, weight1;
    
    if (group0.isEmpty())
    {
        BOOST_ASSERT( not group1.isEmpty() );
        weight0 = 0;
        weight1 = 1;
    }
    else if (group1.isEmpty())
    {
        weight0 = 1;
        weight1 = 0;
    }
    else
    {
        //get the weighting function that is used to weight the two
        //sides of the move
        WeightFunction weightfunc = moldata.property(map["weight function"],
                                                     default_weightfunc);
    
        tuple<double,double> weights = weightfunc(moldata, group0, 
                                                  group1, map);
    
        weight0 = weights.get<0>();
        weight1 = weights.get<1>();
    }
    
    //now get property containing the coordinates of the atoms
    PropertyName coord_property = map["coordinates"];
    
    AtomicCoords coords = moldata.property(coord_property);
    
    //use these coordinates to calculate the unit vector that 
    //points along the bond
    Vector unit_vec = (coords[moldata.info().cgAtomIdx(atom1)] - 
                       coords[moldata.info().cgAtomIdx(atom0)]).normalise();

    //now translate the groups along this vector by their weighted
    //amount of delta
    if (weight0 != 0)
        MoverBase::translate(coords, group0, -weight0 * delta);
        
    if (weight1 != 0)
        MoverBase::translate(coords, group1, weight1 * delta);
        
    //save the new coordinates
    moldata.setProperty(coord_property, coords);
}

/** Change the size of the angle identified by 'angle' by 'delta',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties 
    
    This only moves the movable atoms in this view, and 
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.
    
    The angle is labelled;
    
    atom0     atom2
         \    /
         atom1
    
    The molecule is split by the atom0-..-atom2 bond, i.e. atom0
    is in group0, atom2 is in group1 and atom1 is not in any
    group and is not moved.
    
    The two groups are rotated around the vector perpendicular to
    atom0->atom1 and atom2->atom1, about the point atom1
    
    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::change(MoleculeData &moldata, const AngleID &angle,
                       SireUnits::Dimension::Angle delta,
                       const PropertyMap &map) const
{
    if (delta == 0)
        return;
        
    //get the indicies of the atoms in the angle
    tuple<AtomIdx,AtomIdx,AtomIdx> atomidxs = angle.map(moldata.info());

    AtomIdx atom0 = atomidxs.get<0>();
    AtomIdx atom1 = atomidxs.get<1>();
    AtomIdx atom2 = atomidxs.get<2>();
    
    //get the connectivity that is used to split the 
    //molecule into two parts
    Connectivity connectivity = moldata.property(map["connectivity"]);
    
    //see if there are any anchors that hold part of the 
    //molecule stationary
    AtomSelection anchors;
    
    if (map.specified("anchors"))
        anchors = moldata.property(map["anchors"]);
    
    //split the molecule into the two moving parts
    tuple<AtomSelection,AtomSelection> groups =
                          connectivity.split(atom0, atom1, atom2,
                                             movable_atoms, anchors);
    
    const AtomSelection &group0 = groups.get<0>();
    const AtomSelection &group1 = groups.get<1>();
    
    double weight0, weight1;
    
    if (group0.isEmpty())
    {
        BOOST_ASSERT( not group1.isEmpty() );
        weight0 = 0;
        weight1 = 1;
    }
    else if (group1.isEmpty())
    {
        weight0 = 1;
        weight1 = 0;
    }
    else
    {
        //get the weighting function that is used to weight the 
        //two sides of the move
        WeightFunction weightfunc = moldata.property(map["weight function"],
                                                     default_weightfunc);
        
        tuple<double,double> weights = weightfunc(moldata, group0, 
                                                  group1, map);

        weight0 = weights.get<0>();
        weight1 = weights.get<1>();
    }
    
    //get the coordinates that are to be changed
    PropertyName coord_property = map["coordinates"];
    AtomicCoords coords = moldata.property(coord_property);
    
    //get the coordinates of the three atoms that comprise the angle
    const Vector &coords0 = coords[moldata.info().cgAtomIdx(atom0)];
    const Vector &coords1 = coords[moldata.info().cgAtomIdx(atom1)];
    const Vector &coords2 = coords[moldata.info().cgAtomIdx(atom2)];
    
    //get the vector perpendicular to the angle
    Vector perp = Vector::cross( coords2-coords0, coords1-coords0 );
    
    //rotate the two groups
    if (weight0 != 0)
        MoverBase::rotate(coords, group0, 
                          Quaternion(-weight0*delta, perp), coords1);
                          
    if (weight1 != 0)
        MoverBase::rotate(coords, group1,
                          Quaternion(weight1*delta, perp), coords1);

    //save the new coordinates
    moldata.setProperty(coord_property, coords);
}

/** Change the size of the dihedral identified by 'dihedra' by 'delta',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties 
    
    This only moves the movable atoms in this view, and 
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.
    
    The dihedral is labelled;
    
    atom0           atom3
        \           /
        atom1--atom2
    
    The molecule is split by the atom0-..-atom3 bond (i.e.
    atom0 is in group0, atom3 is in group1, while atom1 and atom2
    are not in any group and are not moved).
    
    The groups are then rotated about the vector atom1->atom2
    
    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::change(MoleculeData &moldata, const DihedralID &dihedral,
                       SireUnits::Dimension::Length delta,
                       const PropertyMap &map) const
{
    if (delta == 0)
        return;
        
    //get the indicies of the atoms that comprise this dihedral
    tuple<AtomIdx,AtomIdx,AtomIdx,AtomIdx> atomidxs = 
                                               dihedral.map(moldata.info());
                                               
    AtomIdx atom0 = atomidxs.get<0>();
    AtomIdx atom1 = atomidxs.get<1>();
    AtomIdx atom2 = atomidxs.get<2>();
    AtomIdx atom3 = atomidxs.get<3>();
    
    //now get the connectivity of the molecule
    Connectivity connectivity = moldata.property(map["connectivity"]);
    
    //see if there are any anchor atoms specified
    AtomSelection anchors;
    
    if (map.specified("anchors"))
        anchors = moldata.property(map["anchors"]);
        
    //split the molecule into the two moving parts
    tuple<AtomSelection,AtomSelection> groups = 
                          connectivity.split(moldata, atom0, atom1,
                                             atom2, atom3,
                                             movable_atoms, anchors);
                                             
    const AtomSelection &group0 = groups.get<0>();
    const AtomSelection &group1 = groups.get<1>();
    
    //get the weights applied to the motion of each group
    double weight0, weight1;
    
    if (group0.isEmpty())
    {
        BOOST_ASSERT( not group1.isEmpty() );
        weight0 = 0;
        weight1 = 1;
    }
    else if (group1.isEmpty())
    {
        weight0 = 1;
        weight1 = 0;
    }
    else
    {
        WeightFunction weightfunc = moldata.property(map["weight function"],
                                                     default_weightfunc);
        
        tuple<double,double> weights = weightfunc(moldata, group0,
                                                  group1, map);
                                                  
        weight0 = weights.get<0>();
        weight1 = weights.get<1>();
    }
    
    //get the coordinates to be moved
    PropertyName coord_property = map["coordinates"];
    AtomicCoords coords = moldata.property(coord_property);
    
    //get the coordinates of the central two atoms of the dihedral
    const Vector &coords1 = coords[moldata.info().cgAtomIdx(atom1)];
    const Vector &coords2 = coords[moldata.info().cgAtomIdx(atom2)];
    
    //get the vector about which the two parts of the molecule
    //are rotated
    Vector dihvec = coords2 - coords1;
    
    //now rotate the two parts of the molecule
    if (weight0 != 0)
        MoverBase::rotate(coords, group0, 
                          Quaternion(-weight*delta, dihvec), coords1);
                          
    if (weight1 != 0)
        MoverBase::rotate(coords, group1,
                          Quaternion(weight*delta, dihvec), coords2);
                          
    //save the new coordinates
    moldata.setProperty(coord_property, coords);
}

/** Change the size of the dihedral identified by 'bond' by 'delta',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties 
    
    This only moves the movable atoms in this view, and 
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.
    
    The dihedral is labelled;
    
      X--atom0--atom1--X
      
    The molecule is split by the atom0-atom1 bond, and the two
    groups are rotated about the vector connecting atom0->atom1
    about the point atom0 for group0, and atom1 for group1
    
    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::change(MoleculeData &moldata, const BondID &bond,
                       SireUnits::Dimension::Angle delta,
                       const PropertyMap &map) const
{
    if (delta == 0)
        return;
        
    //get the indicies of the atoms that comprise this dihedral
    tuple<AtomIdx,AtomIdx> atomidxs = bond.map(moldata.info());
                                               
    AtomIdx atom0 = atomidxs.get<0>();
    AtomIdx atom1 = atomidxs.get<1>();
    
    //now get the connectivity of the molecule
    Connectivity connectivity = moldata.property(map["connectivity"]);
    
    //see if there are any anchor atoms specified
    AtomSelection anchors;
    
    if (map.specified("anchors"))
        anchors = moldata.property(map["anchors"]);
        
    //split the molecule into the two moving parts
    tuple<AtomSelection,AtomSelection> groups = 
                          connectivity.split(moldata, atom0, atom1,
                                             movable_atoms, anchors);
                                             
    const AtomSelection &group0 = groups.get<0>();
    const AtomSelection &group1 = groups.get<1>();
    
    //get the weights applied to the motion of each group
    double weight0, weight1;
    
    if (group0.isEmpty())
    {
        BOOST_ASSERT( not group1.isEmpty() );
        weight0 = 0;
        weight1 = 1;
    }
    else if (group1.isEmpty())
    {
        weight0 = 1;
        weight1 = 0;
    }
    else
    {
        WeightFunction weightfunc = moldata.property(map["weight function"],
                                                     default_weightfunc);
        
        tuple<double,double> weights = weightfunc(moldata, group0,
                                                  group1, map);
                                                  
        weight0 = weights.get<0>();
        weight1 = weights.get<1>();
    }
    
    //get the coordinates to be moved
    PropertyName coord_property = map["coordinates"];
    AtomicCoords coords = moldata.property(coord_property);
    
    //get the coordinates of the central two atoms of the dihedral
    const Vector &coords0 = coords[moldata.info().cgAtomIdx(atom0)];
    const Vector &coords1 = coords[moldata.info().cgAtomIdx(atom1)];
    
    //get the vector about which the two parts of the molecule
    //are rotated
    Vector dihvec = coords1 - coords0;
    
    //now rotate the two parts of the molecule
    if (weight0 != 0)
        MoverBase::rotate(coords, group0, 
                          Quaternion(-weight*delta, dihvec), coords0);
                          
    if (weight1 != 0)
        MoverBase::rotate(coords, group1,
                          Quaternion(weight*delta, dihvec), coords1);
                          
    //save the new coordinates
    moldata.setProperty(coord_property, coords);
}

/** Change the size of the improper identified by 'improper' by 'delta',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties 
    
    This only moves the movable atoms in this view, and 
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.
    
    An improper is labelled;
    
                 atom2
                /
    atom0--atom1
                \
                 atom3
                 
    The molecule is split into two along the atom0-atom1 bond
    (i.e. atom0 is in one group, while atom1, atom2 and atom3 are
    in the other). The groups are then rotated around the vector
    from atom2->atom3, about the point 'atom1'
    
    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::change(MoleculeData &moldata, const ImproperID &improper,
                       SireUnits::Dimension::Angle delta,
                       const PropertyMap &map) const
{
    if (delta == 0)
        return;
        
    //get the indicies of the atoms that comprise this dihedral
    tuple<AtomIdx,AtomIdx,AtomIdx,AtomIdx> atomidxs = 
                                               improper.map(moldata.info());
                                               
    AtomIdx atom0 = atomidxs.get<0>();
    AtomIdx atom1 = atomidxs.get<1>();
    AtomIdx atom2 = atomidxs.get<2>();
    AtomIdx atom3 = atomidxs.get<3>();
    
    //now get the connectivity of the molecule
    Connectivity connectivity = moldata.property(map["connectivity"]);
    
    //see if there are any anchor atoms specified
    AtomSelection anchors;
    
    if (map.specified("anchors"))
        anchors = moldata.property(map["anchors"]);
        
    //split the molecule into the two moving parts
    tuple<AtomSelection,AtomSelection> groups = 
                          connectivity.split(moldata, atom0, atom1,
                                             movable_atoms, anchors);
                                             
    const AtomSelection &group0 = groups.get<0>();
    const AtomSelection &group1 = groups.get<1>();
    
    //get the weights applied to the motion of each group
    double weight0, weight1;
    
    if (group0.isEmpty())
    {
        BOOST_ASSERT( not group1.isEmpty() );
        weight0 = 0;
        weight1 = 1;
    }
    else if (group1.isEmpty())
    {
        weight0 = 1;
        weight1 = 0;
    }
    else
    {
        WeightFunction weightfunc = moldata.property(map["weight function"],
                                                     default_weightfunc);
        
        tuple<double,double> weights = weightfunc(moldata, group0,
                                                  group1, map);
                                                  
        weight0 = weights.get<0>();
        weight1 = weights.get<1>();
    }
    
    //get the coordinates to be moved
    PropertyName coord_property = map["coordinates"];
    AtomicCoords coords = moldata.property(coord_property);
    
    //get the coordinates of the last three atoms of the improper
    const Vector &coords1 = coords[moldata.info().cgAtomIdx(atom1)];
    const Vector &coords2 = coords[moldata.info().cgAtomIdx(atom2)];
    const Vector &coords3 = coords[moldata.info().cgAtomIdx(atom3)];
    
    //get the vector from atom2 to atom3
    Vector impvec = coords3 - coords2;
    
    //now rotate the two parts of the molecule
    if (weight0 != 0)
        MoverBase::rotate(coords, group0, 
                          Quaternion(-weight*delta, impvec), coords1);
                          
    if (weight1 != 0)
        MoverBase::rotate(coords, group1,
                          Quaternion(weight*delta, impvec), coords1);
                          
    //save the new coordinates
    moldata.setProperty(coord_property, coords);
}

/** Set the length of the bond identified by 'bond' to 'value',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties 
    
    This only moves the movable atoms in this view, and 
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.
    
    The bond is labelled;
    
    atom0--atom1
    
    The molecule is split into two about this bond, i.e.
    atom0 and everything it is bonded to is in group0, while
    atom1 and everything it is bonded to is in group1.
    
    The two groups are then translated along the vector atom0->atom1
    
    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::set(MoleculeData &moldata, const BondID &bond,
                    SireUnits::Dimension::Length value,
                    const PropertyMap &map) const
{
    double current_value = bond.size(moldata,map);
    
    this->change(moldata, bond, value - current_value, map);
}

/** Set the size of the angle identified by 'angle' to 'value',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties 
    
    This only moves the movable atoms in this view, and 
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.
    
    The angle is labelled;
    
    atom0     atom2
         \    /
         atom1
    
    The molecule is split by the atom0-..-atom2 bond, i.e. atom0
    is in group0, atom2 is in group1 and atom1 is not in any
    group and is not moved.
    
    The two groups are rotated around the vector perpendicular to
    atom0->atom1 and atom2->atom1, about the point atom1
    
    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::set(MoleculeData &moldata, const AngleID &angle,
                    SireUnits::Dimension::Angle value,
                    const PropertyMap &map) const
{
    SireUnits::Dimension::Angle current_value = angle.size(moldata, map);
    
    this->change(moldata, angle, current_value - value, map);
}

/** Set the size of the dihedral identified by 'dihedra' to 'value',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties 
    
    This only moves the movable atoms in this view, and 
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.
    
    The dihedral is labelled;
    
    atom0           atom3
        \           /
        atom1--atom2
    
    The molecule is split by the atom0-..-atom3 bond (i.e.
    atom0 is in group0, atom3 is in group1, while atom1 and atom2
    are not in any group and are not moved).
    
    The groups are then rotated about the vector atom1->atom2
    
    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::set(MoleculeData &moldata, const DihedralID &dihedral,
                    SireUnits::Dimension::Length value,
                    const PropertyMap &map) const
{
    SireUnits::Dimension::Angle current_value = dihedral.size(moldata, map);
    
    this->change(moldata, dihedral, value - current_value, map);
}

/** Set the size of the dihedral identified by 'dihedral' to 'value',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties 
    
    This only moves the movable atoms in this view, and 
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.
    
    The dihedral is labelled;
    
    atom0           atom3
        \           /
        atom1--atom2
    
    The molecule is split by the atom0-..-atom3 bond (i.e.
    atom0 is in group0, atom3 is in group1, while atom1 and atom2
    are not in any group and are not moved).
    
    The groups are then rotated about the vector atom1->atom2
    
    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::setAll(MoleculeData &moldata, const DihedralID &dihedral,
                       SireUnits::Dimension::Angle value,
                       const PropertyMap &map) const
{
    SireUnits::Dimension::Angle current_value = dihedral.size(moldata, map);
    
    this->change(moldata, dihedral.bond12(), 
                 value - current_value, map);
}

/** Set the size of the improper identified by 'improper' to 'value',
    in the molecule whose data is in 'moldata', using the supplied
    PropertyMap to locate the necessary properties 
    
    This only moves the movable atoms in this view, and 
    an anchor_error is thrown if it is not possible to make
    this change without moving the unmovable atoms.
    
    An improper is labelled;
    
                 atom2
                /
    atom0--atom1
                \
                 atom3
                 
    The molecule is split into two along the atom0-atom1 bond
    (i.e. atom0 is in one group, while atom1, atom2 and atom3 are
    in the other). The groups are then rotated around the vector
    from atom2->atom3, about the point 'atom1'
    
    \throw SireBase::missing_property
    \throw SireMol::anchor_error
    \throw SireMol::ring_error
*/
void MoverBase::set(MoleculeData &moldata, const ImproperID &improper,
                    SireUnits::Dimension::Angle value,
                    const PropertyMap &map) const
{
    SireUnits::Dimension::Angle current_value = improper.size(moldata, map);
    
    this->change(moldata, improper, value - current_value, map);
}

/** Align the movable atoms of 'data' against their equivalents in 'other',
    using the 'aliases' where necessary to find the equivalent atom.
    
    If this can't find an equivalent atom, then that atom is ignored.
    If no atoms are found at all, then an exception is raised
    
    This only moves the movable atoms in this view.
    
    This works by calculating the alignment axes and then
    mapping the movable atoms from their current frame into 
    the alignment axes frame
    
    \throw SireMol::missing_atom
*/  
void MoverBase::align(MoleculeData &moldata, const MoleculeView &other,
                      const AtomAliases &aliases,
                      const PropertyMap &map0, const PropertyMap &map1) const
{
    AxisSet alignment_axes = Evaluator::alignmentAxes(moldata, movable_atoms,
                                                      other, aliases,
                                                      map0, map1);
                                                      
    MoverBase::mapInto(moldata, movable_atoms, alignment_axes, map0);
}

/** Align the movable atoms of 'data', using the alignment axes
    calculated to align 'align_atoms' of this molecule with their
    equivalens in 'other', using 'aliases' if necessary to map
    from the atom ID in this molecule to the atom ID in 'other'.
    
    If an equivalent atom can't be found then it is ignored - if
    not equivalent atoms are found then an exception is raised.
    
    This only moves the movable atoms in this view.
    
    This works by calculating the alignment axes and then
    mapping the movable atoms from their current frame into
    the alignment axes frame
    
    \throw SireMol::missing_atom
*/  
void MoverBase::align(MoleculeData &data, const AtomSelection &align_atoms,
                      const MoleculeView &other, const AtomAliases &aliases,
                      const PropertyMap &map0, const PropertyMap &map1) const
{
    AxisSet alignment_axes = Evaluator::alignmentAxes(moldata, align_atoms,
                                                      other, aliases,
                                                      map0, map1);

    MoverBase::mapInto(moldata, movable_atoms, alignment_axes, map0);
}
