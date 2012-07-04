// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License


#include "boost/python.hpp"

#include "Helpers/clone_const_reference.hpp"

#include "DLMRigidBody.pypp.hpp"

#include "DofID.pypp.hpp"

#include "Dynamics.pypp.hpp"

#include "Ensemble.pypp.hpp"

#include "Flexibility.pypp.hpp"

#include "GetCOGPoint.pypp.hpp"

#include "GetCOMPoint.pypp.hpp"

#include "GetPoint.pypp.hpp"

#include "HMCGenerator.pypp.hpp"

#include "HMCVelGen.pypp.hpp"

#include "HybridMC.pypp.hpp"

#include "Integrator.pypp.hpp"

#include "InternalMove.pypp.hpp"

#include "InternalMoveSingle.pypp.hpp"

#include "MTSMC.pypp.hpp"

#include "MaxwellBoltzmann.pypp.hpp"

#include "MolDeleter.pypp.hpp"

#include "MolInserter.pypp.hpp"

#include "MolecularDynamics.pypp.hpp"

#include "MonteCarlo.pypp.hpp"

#include "Move.pypp.hpp"

#include "Moves.pypp.hpp"

#include "NullDeleter.pypp.hpp"

#include "NullGetPoint.pypp.hpp"

#include "NullInserter.pypp.hpp"

#include "NullIntegrator.pypp.hpp"

#include "NullMove.pypp.hpp"

#include "NullSupraMove.pypp.hpp"

#include "NullSupraSubMove.pypp.hpp"

#include "NullVelocityGenerator.pypp.hpp"

#include "NullVolumeChanger.pypp.hpp"

#include "OpenMMIntegrator.pypp.hpp"

#include "PrefSampler.pypp.hpp"

#include "RBWorkspaceJM.pypp.hpp"

#include "RepExMove.pypp.hpp"

#include "RepExSubMove.pypp.hpp"

#include "Replica.pypp.hpp"

#include "Replicas.pypp.hpp"

#include "RigidBodyMC.pypp.hpp"

#include "SameMoves.pypp.hpp"

#include "SameSupraMoves.pypp.hpp"

#include "SameSupraSubMoves.pypp.hpp"

#include "Sampler.pypp.hpp"

#include "ScaleVolumeFromCenter.pypp.hpp"

#include "SimPacket.pypp.hpp"

#include "SimStore.pypp.hpp"

#include "Simulation.pypp.hpp"

#include "SpecifiedGroupsDeleter.pypp.hpp"

#include "SupraMove.pypp.hpp"

#include "SupraMoves.pypp.hpp"

#include "SupraSim.pypp.hpp"

#include "SupraSimPacket.pypp.hpp"

#include "SupraSubMove.pypp.hpp"

#include "SupraSubMoves.pypp.hpp"

#include "SupraSubSim.pypp.hpp"

#include "SupraSubSimPacket.pypp.hpp"

#include "SupraSubSystem.pypp.hpp"

#include "SupraSystem.pypp.hpp"

#include "SystemWideDeleter.pypp.hpp"

#include "UniformInserter.pypp.hpp"

#include "UniformSampler.pypp.hpp"

#include "VelocitiesFromProperty.pypp.hpp"

#include "VelocityGenerator.pypp.hpp"

#include "VelocityVerlet.pypp.hpp"

#include "VelocityVerletBerendsen.pypp.hpp"

#include "VolumeChanger.pypp.hpp"

#include "VolumeMove.pypp.hpp"

#include "WeightedMoves.pypp.hpp"

#include "ZMatMove.pypp.hpp"

#include "ZMatrix.pypp.hpp"

#include "ZMatrixCoords.pypp.hpp"

#include "ZMatrixCoordsLine.pypp.hpp"

#include "ZMatrixLine.pypp.hpp"

namespace bp = boost::python;

#include "SireMove_containers.h"

#include "SireMove_registrars.h"

#include "SireMove_properties.h"

BOOST_PYTHON_MODULE(_Move){
    register_SireMove_objects();

    register_SireMove_containers();

    register_Integrator_class();

    register_DLMRigidBody_class();

    register_DofID_class();

    register_Move_class();

    register_Dynamics_class();

    register_Ensemble_class();

    register_Flexibility_class();

    register_GetPoint_class();

    register_GetCOGPoint_class();

    register_GetCOMPoint_class();

    register_VelocityGenerator_class();

    register_HMCVelGen_class();

    register_HMCGenerator_class();

    register_MonteCarlo_class();

    register_HybridMC_class();

    register_InternalMove_class();

    register_InternalMoveSingle_class();

    register_MTSMC_class();

    register_MaxwellBoltzmann_class();

    register_MolDeleter_class();

    register_MolInserter_class();

    register_MolecularDynamics_class();

    register_Moves_class();

    register_NullDeleter_class();

    register_NullGetPoint_class();

    register_NullInserter_class();

    register_NullIntegrator_class();

    register_NullMove_class();

    register_SupraMove_class();

    register_NullSupraMove_class();

    register_SupraSubMove_class();

    register_NullSupraSubMove_class();

    register_NullVelocityGenerator_class();

    register_VolumeChanger_class();

    register_NullVolumeChanger_class();

    register_OpenMMIntegrator_class();

    register_Sampler_class();

    register_PrefSampler_class();

    register_RBWorkspaceJM_class();

    register_RepExMove_class();

    register_RepExSubMove_class();

    register_SupraSubSystem_class();

    register_Replica_class();

    register_SupraSystem_class();

    register_Replicas_class();

    register_RigidBodyMC_class();

    register_SameMoves_class();

    register_SupraMoves_class();

    register_SameSupraMoves_class();

    register_SupraSubMoves_class();

    register_SameSupraSubMoves_class();

    register_ScaleVolumeFromCenter_class();

    register_SimPacket_class();

    register_SimStore_class();

    register_Simulation_class();

    register_SpecifiedGroupsDeleter_class();

    register_SupraSim_class();

    register_SupraSimPacket_class();

    register_SupraSubSim_class();

    register_SupraSubSimPacket_class();

    register_SystemWideDeleter_class();

    register_UniformInserter_class();

    register_UniformSampler_class();

    register_VelocitiesFromProperty_class();

    register_VelocityVerlet_class();

    register_VelocityVerletBerendsen_class();

    register_VolumeMove_class();

    register_WeightedMoves_class();

    register_ZMatMove_class();

    register_ZMatrix_class();

    register_ZMatrixCoords_class();

    register_ZMatrixLine_class();

    register_SireMove_properties();

    register_ZMatrixCoordsLine_class();
}

