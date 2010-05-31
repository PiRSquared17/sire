/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#include "dlmrigidbody.h"
#include "ensemble.h"
#include "rbworkspace.h"

#include "SireMaths/quaternion.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"

#include "SireStream/datastream.h"

using namespace SireMove;
using namespace SireMol;
using namespace SireSystem;
using namespace SireFF;
using namespace SireMaths;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<DLMRigidBody> r_dlmrb;

QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const DLMRigidBody &dlmrb)
{
    writeHeader(ds, r_dlmrb, 1);
    
    ds << dlmrb.frequent_save_velocities << static_cast<const Integrator&>(dlmrb);
    
    return ds;
}

QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, DLMRigidBody &dlmrb)
{
    VersionID v = readHeader(ds, r_dlmrb);
    
    if (v == 1)
    {
        ds >> dlmrb.frequent_save_velocities >> static_cast<Integrator&>(dlmrb);
    }
    else
        throw version_error( v, "1", r_dlmrb, CODELOC );
        
    return ds;
}

/** Constructor */
DLMRigidBody::DLMRigidBody(bool frequent_save) 
             : ConcreteProperty<DLMRigidBody,Integrator>(),
               frequent_save_velocities(frequent_save)
{}

/** Copy constructor */
DLMRigidBody::DLMRigidBody(const DLMRigidBody &other)
             : ConcreteProperty<DLMRigidBody,Integrator>(other),
               frequent_save_velocities(other.frequent_save_velocities)
{}

/** Destructor */
DLMRigidBody::~DLMRigidBody()
{}

/** Copy assignment operator */
DLMRigidBody& DLMRigidBody::operator=(const DLMRigidBody &other)
{
    if (this != &other)
    {
        Integrator::operator=(other);
        frequent_save_velocities = other.frequent_save_velocities;
    }
        
    return *this;
}

/** Comparison operator */
bool DLMRigidBody::operator==(const DLMRigidBody &other) const
{
    return frequent_save_velocities == other.frequent_save_velocities and 
           Integrator::operator==(other);
}

/** Comparison operator */
bool DLMRigidBody::operator!=(const DLMRigidBody &other) const
{
    return not DLMRigidBody::operator==(other);
}

const char* DLMRigidBody::typeName()
{
    return QMetaType::typeName( qMetaTypeId<DLMRigidBody>() );
}

QString DLMRigidBody::toString() const
{
    return QObject::tr("DLMRigidBody()");
}

/** Return the ensemble generated by this integrator (NVE) */
Ensemble DLMRigidBody::ensemble() const
{
    return Ensemble::NVE();
}

/** This integrator is time-reversible */
bool DLMRigidBody::isTimeReversible() const
{
    return true;
}

static Matrix getRX(double phi)
{
    return Quaternion( phi*radian, Vector(1,0,0) ).toMatrix();

    return Matrix( 1, 0, 0,
                   0, cos(phi), -sin(phi),
                   0, sin(phi), cos(phi) );
}

static Matrix getRY(double phi)
{
    return Quaternion( phi*radian, Vector(0,1,0) ).toMatrix();

    return Matrix( cos(phi), 0, sin(phi),
                   0, 1, 0,
                   -sin(phi), 0, cos(phi) );
}

static Matrix getRZ(double phi)
{
    return Quaternion( phi*radian, Vector(0,0,1) ).toMatrix();

    return Matrix( cos(phi), -sin(phi), 0,
                   sin(phi), cos(phi), 0,
                   0, 0, 1 );
}

/** Integrate using the passed workspace, the energy component 'nrg_component',
    performing 'nmoves' moves using the specified timestep, recording
    statistics of every move if 'record_stats' is true */
void DLMRigidBody::integrate(IntegratorWorkspace &workspace,
                             const Symbol &nrg_component, 
                             Time timestep, int nmoves, bool record_stats) const
{
    RBWorkspace &ws = workspace.asA<RBWorkspace>();
    
    const double dt = timestep.value();
    const double half_dt = 0.5 * dt;

    const int nbeads = ws.nBeads();
    
    const Vector X(1,0,0);
    const Vector Y(0,1,0);
    const Vector Z(0,0,1);
    
    for (int imove=0; imove<nmoves; ++imove)
    {
        ws.calculateForces(nrg_component);
        
        Vector *bead_coords = ws.beadCoordsArray();
        Quaternion *bead_orient = ws.beadOrientationArray();
        Vector *bead_lin_momenta = ws.beadLinearMomentaArray();
        Vector *bead_ang_momenta = ws.beadAngularMomentaArray();
        const double *bead_masses = ws.beadMassesArray();
        const Vector *bead_inertia = ws.beadInertiasArray();
        
        const Vector *bead_forces = ws.beadForcesArray();
        const Vector *bead_torques = ws.beadTorquesArray();
        
        //first integrate the coordinates - loop over all beads
        for (int i=0; i<nbeads; ++i)
        {
            Vector &x = bead_coords[i];
            Quaternion &q = bead_orient[i];
            Vector &p = bead_lin_momenta[i];
            Vector &ap = bead_ang_momenta[i];
            double mass = bead_masses[i];
            const Vector &inertia = bead_inertia[i];
            const Vector &force = bead_forces[i];
            const Vector &torque = bead_torques[i];

            if (mass == 0)
                //this is a dummy atom
                continue;

            // use velocity verlet to integrate the position of the bead
            // v(t + dt/2) = v(t) + (1/2) a(t) dt
            p += (half_dt * force);

            // r(t + dt) = r(t) + v(t + dt/2) dt
            x += (dt / mass) * p;
            
            //now update the orientation / angular momenta using the DLM algorithm
            ap += (half_dt * torque);
            
            if (not ap.isZero())
            {
                if (ap[0] != 0 and inertia[0] != 0)
                {
                    //Matrix R1 = Quaternion( (half_dt * ap[0] / inertia[0])*radian, X )
                    //                .toMatrix();
                    //Matrix R1 = ::getRX(half_dt * ap[0] / inertia[0]);
                    Quaternion R1( (half_dt*ap[0] / inertia[0]) * radian, Vector(1,0,0) );
            
                    ap = R1.rotate(ap);
                    q = q * R1.conjugate();
                }
            
                if (ap[1] != 0 and inertia[1] != 0)
                {
                    //Matrix R2 = Quaternion( (half_dt * ap[1] / inertia[1])*radian, Y )
                    //                 .toMatrix();
                    //Matrix R2 = ::getRY(half_dt * ap[1] / inertia[1]);
                    Quaternion R2( (half_dt*ap[1] / inertia[1]) * radian, Vector(0,1,0) );
                                          
                    ap = R2.rotate(ap);
                    q = q * R2.conjugate();
                }
            
                if (ap[2] != 0 and inertia[2] != 0)
                {
                    //Matrix R3 = Quaternion( (dt * ap[2] / inertia[2])*radian, Z )
                    //                    .toMatrix();
                    //Matrix R3 = ::getRZ(dt * ap[2] / inertia[2]);
                    Quaternion R3( (dt*ap[2] / inertia[2]) * radian, Vector(0,0,1) );
            
                    ap = R3.rotate(ap);
                    q = q * R3.conjugate();
                }
                
                if (ap[1] != 0 and inertia[1] != 0)
                {
                    //Matrix R4 = Quaternion( (half_dt * ap[1] / inertia[1])*radian, Y )
                    //                    .toMatrix();
                    //Matrix R4 = ::getRY(half_dt * ap[1] / inertia[1]);
                    Quaternion R4( (half_dt*ap[1] / inertia[1]) * radian, Vector(0,1,0) );
            
                    ap = R4.rotate(ap);
                    q = q * R4.conjugate();
                }
            
                if (ap[0] != 0 and inertia[0] != 0)
                {
                    //Matrix R5 = Quaternion( (half_dt * ap[0] / inertia[0])*radian, X )
                    //                    .toMatrix();
                    //Matrix R5 = ::getRX(half_dt * ap[0] / inertia[0]);
                    Quaternion R5( (half_dt*ap[0] / inertia[0]) * radian, Vector(1,0,0) );
                                  
                    ap = R5.rotate(ap);
                    q = q * R5.conjugate();
                }
            }
        }

        ws.commitCoordinates();
        ws.calculateForces(nrg_component);

        bead_lin_momenta = ws.beadLinearMomentaArray();
        bead_ang_momenta = ws.beadAngularMomentaArray();
        bead_masses = ws.beadMassesArray();
        bead_inertia = ws.beadInertiasArray();
        
        bead_forces = ws.beadForcesArray();
        bead_torques = ws.beadTorquesArray();
        
        //now need to integrate the momenta
        for (int i=0; i<nbeads; ++i)
        {
            if (bead_masses[i] == 0)
                continue;
        
            Vector &p = bead_lin_momenta[i];
            Vector &ap = bead_ang_momenta[i];
            const Vector &force = bead_forces[i];
            const Vector &torque = bead_torques[i];

            p += ( half_dt * force );
            ap += ( half_dt * torque );
        }
        
        if (frequent_save_velocities)
            ws.commitVelocities();
        
        if (record_stats)
            ws.collectStatistics();
    }
    
    if (not frequent_save_velocities)
        ws.commitVelocities();
}

/** Create an empty workspace for this integrator */
IntegratorWorkspacePtr DLMRigidBody::createWorkspace(const PropertyMap &map) const
{
    return IntegratorWorkspacePtr(new RBWorkspace(map));
}

/** Create a workspace to integrate the molecules in the passed molecule group */
IntegratorWorkspacePtr DLMRigidBody::createWorkspace(const MoleculeGroup &molgroup,
                                                     const PropertyMap &map) const
{
    return IntegratorWorkspacePtr(new RBWorkspace(molgroup,map));
}
