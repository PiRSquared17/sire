########################################
#
# CMake file for library:SireMove
#
########################################

# Third Party dependencies of this module

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

# Define the headers in SireMove
set ( SIREMOVE_HEADERS
      dlmrigidbody.h
      dynamics.h
      ensemble.h
      errors.h
      getpoint.h
      flexibility.h
      gibbsmove.h
      hybridmc.h
      integrator.h
      integratorworkspace.h
      internalmove.h
      internalmovesingle.h
      moldeleter.h
      molinserter.h
      moleculardynamics.h
      montecarlo.h
      move.h
      moves.h
      mtsmc.h
      prefsampler.h
      rbworkspace.h
      repexmove.h
      replica.h
      replicas.h
      rigidbodymc.h
      sampler.h
      simpacket.h
      simstore.h
      simulation.h
      supramove.h
      supramoves.h
      suprasim.h
      suprasimpacket.h
      suprasubmove.h
      suprasubmoves.h
      suprasubsim.h
      suprasubsimpacket.h
      suprasubsystem.h
      suprasystem.h
      titrationmove.h
      titrator.h
      uniformsampler.h
      velocitygenerator.h
      velocityverlet.h
      volumechanger.h
      volumemove.h
      weightedmoves.h
      zmatmove.h
      zmatrix.h
    )

# Define the sources in SireMove
set ( SIREMOVE_SOURCES

      register_siremove.cpp

      dlmrigidbody.cpp
      dynamics.cpp
      ensemble.cpp
      errors.cpp
      getpoint.cpp
      flexibility.cpp
      gibbsmove.cpp
      hybridmc.cpp
      integrator.cpp
      integratorworkspace.cpp
      internalmove.cpp
      internalmovesingle.cpp
      moldeleter.cpp
      molinserter.cpp
      moleculardynamics.cpp
      montecarlo.cpp
      move.cpp
      moves.cpp
      mtsmc.cpp
      prefsampler.cpp
      rbworkspace.cpp
      repexmove.cpp
      replica.cpp
      replicas.cpp
      rigidbodymc.cpp
      sampler.cpp
      simpacket.cpp
      simstore.cpp
      simulation.cpp
      supramove.cpp
      supramoves.cpp
      suprasim.cpp
      suprasimpacket.cpp
      suprasubmove.cpp
      suprasubmoves.cpp
      suprasubsim.cpp
      suprasubsimpacket.cpp
      suprasubsystem.cpp
      suprasystem.cpp
      titrationmove.cpp
      titrator.cpp
      uniformsampler.cpp
      weightedmoves.cpp
      velocitygenerator.cpp
      velocityverlet.cpp
      volumechanger.cpp
      volumemove.cpp
      zmatmove.cpp
      zmatrix.cpp   

      ${SIREMOVE_HEADERS}
    )

add_library (SireMove ${SIREMOVE_SOURCES})

set_target_properties (SireMove PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${SIRE_VERSION_MAJOR}
                      )

# Link to other Sire libraries
target_link_libraries (SireMove
                       SireSystem
                       SireMaths
                       SireCluster
                       SireStream
                       )

# installation
install( TARGETS SireMove EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIREMOVE_HEADERS}
         DESTINATION ${SIRE_INCLUDES}/SireMove )

install( FILES ${SIREMOVE_DEPRECATED_HEADERS}
         DESTINATION ${SIRE_INCLUDES}/SireMove/deprecated )
