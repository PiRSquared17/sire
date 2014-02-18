/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#ifndef SIREMM_CLJSHIFTFUNCTION_H
#define SIREMM_CLJSHIFTFUNCTION_H

#include "cljfunction.h"

#include <QHash>

namespace SireMM
{
class CLJShiftFunction;
class CLJIntraShiftFunction;
class CLJSoftShiftFunction;
class CLJSoftIntraShiftFunction;
}

QDataStream& operator<<(QDataStream&, const SireMM::CLJShiftFunction&);
QDataStream& operator>>(QDataStream&, SireMM::CLJShiftFunction&);

QDataStream& operator<<(QDataStream&, const SireMM::CLJIntraShiftFunction&);
QDataStream& operator>>(QDataStream&, SireMM::CLJIntraShiftFunction&);

namespace SireMM
{

/** This CLJFunction calculates the intermolecular coulomb and LJ energy of the passed
    CLJAtoms using a force shifted electrostatics cutoff
    
    @author Christopher Woods
*/
class SIREMM_EXPORT CLJShiftFunction
        : public SireBase::ConcreteProperty<CLJShiftFunction,CLJCutoffFunction>
{

friend QDataStream& ::operator<<(QDataStream&, const CLJShiftFunction&);
friend QDataStream& ::operator>>(QDataStream&, CLJShiftFunction&);

public:
    CLJShiftFunction();
    CLJShiftFunction(Length cutoff);
    CLJShiftFunction(Length coul_cutoff, Length lj_cutoff);
    
    CLJShiftFunction(const Space &space, Length cutoff);
    CLJShiftFunction(const Space &space, Length coul_cutoff, Length lj_cutoff);
    
    CLJShiftFunction(Length cutoff, COMBINING_RULES combining_rules);
    CLJShiftFunction(Length coul_cutoff, Length lj_cutoff, COMBINING_RULES combining_rules);
    
    CLJShiftFunction(const Space &space, COMBINING_RULES combining_rules);
    CLJShiftFunction(const Space &space, Length cutoff, COMBINING_RULES combining_rules);
    CLJShiftFunction(const Space &space, Length coul_cutoff, Length lj_cutoff,
                     COMBINING_RULES combining_rules);
    
    CLJShiftFunction(const CLJShiftFunction &other);
    
    ~CLJShiftFunction();
    
    CLJShiftFunction& operator=(const CLJShiftFunction &other);
    
    bool operator==(const CLJShiftFunction &other) const;
    bool operator!=(const CLJShiftFunction &other) const;
    
    static const char* typeName();
    const char* what() const;
    
    CLJShiftFunction* clone() const;

    bool supportsGridCalculation() const;

protected:
    void calcVacEnergyAri(const CLJAtoms &atoms,
                          double &cnrg, double &ljnrg) const;
    
    void calcVacEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          double &cnrg, double &ljnrg, float min_distance) const;

    void calcVacEnergyGeo(const CLJAtoms &atoms,
                          double &cnrg, double &ljnrg) const;
    
    void calcVacEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          double &cnrg, double &ljnrg, float min_distance) const;

    void calcBoxEnergyAri(const CLJAtoms &atoms, const Vector &box_dimensions,
                          double &cnrg, double &ljnrg) const;
    
    void calcBoxEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          const Vector &box_dimensions, double &cnrg, double &ljnrg,
                          float min_distance) const;

    void calcBoxEnergyGeo(const CLJAtoms &atoms, const Vector &box_dimensions,
                          double &cnrg, double &ljnrg) const;
    
    void calcBoxEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          const Vector &box_dimensions, double &cnrg, double &ljnrg,
                          float min_distance) const;

    double calcVacCoulombEnergyAri(const CLJAtoms &atoms) const;
    double calcVacCoulombEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                   float min_distance) const;
    
    double calcVacLJEnergyAri(const CLJAtoms &atoms) const;
    double calcVacLJEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                              float min_distance) const;

    void calcVacGrid(const CLJAtoms &atoms, const GridInfo &gridinfo,
                     const int start, const int end, float *gridpot) const;
    
    void calcBoxGrid(const CLJAtoms &atoms, const GridInfo &gridinfo,
                     const int start, const int end, float *gridpot) const;
};

/** This CLJFunction calculates the intramolecular coulomb and LJ energy of the passed
    CLJAtoms using a force shifted electrostatics cutoff
    
    @author Christopher Woods
*/
class SIREMM_EXPORT CLJIntraShiftFunction
        : public SireBase::ConcreteProperty<CLJIntraShiftFunction,CLJIntraFunction>
{

friend QDataStream& ::operator<<(QDataStream&, const CLJIntraShiftFunction&);
friend QDataStream& ::operator>>(QDataStream&, CLJIntraShiftFunction&);

public:
    CLJIntraShiftFunction();
    CLJIntraShiftFunction(Length cutoff);
    CLJIntraShiftFunction(Length coul_cutoff, Length lj_cutoff);
    
    CLJIntraShiftFunction(const Space &space, Length cutoff);
    CLJIntraShiftFunction(const Space &space, Length coul_cutoff, Length lj_cutoff);
    
    CLJIntraShiftFunction(Length cutoff, COMBINING_RULES combining_rules);
    CLJIntraShiftFunction(Length coul_cutoff, Length lj_cutoff, COMBINING_RULES combining_rules);
    
    CLJIntraShiftFunction(const Space &space, COMBINING_RULES combining_rules);
    CLJIntraShiftFunction(const Space &space, Length cutoff, COMBINING_RULES combining_rules);
    CLJIntraShiftFunction(const Space &space, Length coul_cutoff, Length lj_cutoff,
                          COMBINING_RULES combining_rules);
    
    CLJIntraShiftFunction(const CLJIntraShiftFunction &other);
    
    ~CLJIntraShiftFunction();
    
    CLJIntraShiftFunction& operator=(const CLJIntraShiftFunction &other);
    
    bool operator==(const CLJIntraShiftFunction &other) const;
    bool operator!=(const CLJIntraShiftFunction &other) const;
    
    static const char* typeName();
    const char* what() const;
    
    CLJIntraShiftFunction* clone() const;

protected:
    void calcVacEnergyAri(const CLJAtoms &atoms,
                          double &cnrg, double &ljnrg) const;
    
    void calcVacEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          double &cnrg, double &ljnrg, float min_distance) const;

    void calcVacEnergyGeo(const CLJAtoms &atoms,
                          double &cnrg, double &ljnrg) const;
    
    void calcVacEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          double &cnrg, double &ljnrg, float min_distance) const;

    void calcBoxEnergyAri(const CLJAtoms &atoms, const Vector &box_dimensions,
                          double &cnrg, double &ljnrg) const;
    
    void calcBoxEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          const Vector &box_dimensions, double &cnrg, double &ljnrg
                          , float min_distance) const;

    void calcBoxEnergyGeo(const CLJAtoms &atoms, const Vector &box_dimensions,
                          double &cnrg, double &ljnrg) const;
    
    void calcBoxEnergyGeo(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                          const Vector &box_dimensions, double &cnrg, double &ljnrg,
                          float min_distance) const;

    double calcVacCoulombEnergyAri(const CLJAtoms &atoms) const;
    double calcVacCoulombEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                                   float min_distance) const;
    
    double calcVacLJEnergyAri(const CLJAtoms &atoms) const;
    double calcVacLJEnergyAri(const CLJAtoms &atoms0, const CLJAtoms &atoms1,
                              float min_distance) const;
};

}

Q_DECLARE_METATYPE( SireMM::CLJShiftFunction )
Q_DECLARE_METATYPE( SireMM::CLJIntraShiftFunction )

SIRE_EXPOSE_CLASS( SireMM::CLJShiftFunction )
SIRE_EXPOSE_CLASS( SireMM::CLJIntraShiftFunction )

#endif
