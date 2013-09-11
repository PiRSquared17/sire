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

#ifndef SIREMATHS_FREEENERGYAVERAGE_H
#define SIREMATHS_FREEENERGYAVERAGE_H

#include "accumulator.h"
#include "histogram.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/temperature.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
class FreeEnergyAverage;
class BennettsFreeEnergyAverage;
}

QDataStream& operator<<(QDataStream&, const SireMaths::FreeEnergyAverage&);
QDataStream& operator>>(QDataStream&, SireMaths::FreeEnergyAverage&);

QDataStream& operator<<(QDataStream&, const SireMaths::BennettsFreeEnergyAverage&);
QDataStream& operator>>(QDataStream&, SireMaths::BennettsFreeEnergyAverage&);

namespace SireMaths
{

/** This class is used to accumulate the free energy
    average. As well as calculating the average, it also
    records a histogram of values that can be used for
    error analysis
    
    @author Christopher Woods
*/
class SIREMATHS_EXPORT FreeEnergyAverage
           : public SireBase::ConcreteProperty<FreeEnergyAverage,ExpAverage>
{

friend QDataStream& ::operator<<(QDataStream&, const FreeEnergyAverage&);
friend QDataStream& ::operator>>(QDataStream&, FreeEnergyAverage&);

public:
    FreeEnergyAverage();
    FreeEnergyAverage(const SireUnits::Dimension::Temperature &temperature);
    
    FreeEnergyAverage(const SireUnits::Dimension::MolarEnergy &binwidth);
    
    FreeEnergyAverage(const SireUnits::Dimension::Temperature &temperature,
                      const SireUnits::Dimension::MolarEnergy &binwidth);
    
    FreeEnergyAverage(const FreeEnergyAverage &other);
    
    ~FreeEnergyAverage();
    
    FreeEnergyAverage& operator=(const FreeEnergyAverage &other);
    
    bool operator==(const FreeEnergyAverage &other) const;
    bool operator!=(const FreeEnergyAverage &other) const;
    
    FreeEnergyAverage operator+(const FreeEnergyAverage &other) const;
    
    FreeEnergyAverage& operator+=(const FreeEnergyAverage &other);
    
    static const char* typeName();

    SireUnits::Dimension::Temperature temperature() const;

    QString toString() const;

    const Histogram& histogram() const;

    void clear();
    
    void accumulate(double value);

    double taylorExpansion() const;

private:
    /** A histogram of all of the energy values. This allows the free energy
        average to be monitored and the error, and corrected values to be 
        calculated */
    Histogram hist;
};

/** This class is used to accumulate the free energy average, using
    both FEP and Bennett's acceptance ratio method. Bennett's method
    calculates a free energy difference between A and B using 
    simulations at both ensembles A and B. The forwards average
    < f( beta(U_B-U_A) )>_A divided by the backwards average
    < f( beta(U_A-U_B) )>_B is equal to exp(-beta deltaG).
    
    In this case, we don't use the "C" energy offset value, as this
    is optimally equal to deltaG, which we don't know, and also because
    in windowed calculations, abs(deltaG) is < 5 kcal mol-1.
    
    Also note that we use the function f(x) = 1 / (1 + exp(x)) as
    this is also supposed to be the most efficient. 
 
    @author Christopher Woods
*/
class SIREMATHS_EXPORT BennettsFreeEnergyAverage
       : public SireBase::ConcreteProperty<BennettsFreeEnergyAverage,FreeEnergyAverage>
{

friend QDataStream& ::operator<<(QDataStream&, const BennettsFreeEnergyAverage&);
friend QDataStream& ::operator>>(QDataStream&, BennettsFreeEnergyAverage&);

public:
    BennettsFreeEnergyAverage();

    BennettsFreeEnergyAverage(const SireUnits::Dimension::Temperature &temperature);
    
    BennettsFreeEnergyAverage(const SireUnits::Dimension::MolarEnergy &binwidth);
    
    BennettsFreeEnergyAverage(const SireUnits::Dimension::Temperature &temperature,
                             const SireUnits::Dimension::MolarEnergy &binwidth);
    
    BennettsFreeEnergyAverage(const BennettsFreeEnergyAverage &other);
    
    ~BennettsFreeEnergyAverage();
    
    BennettsFreeEnergyAverage& operator=(const BennettsFreeEnergyAverage &other);
    
    bool operator==(const BennettsFreeEnergyAverage &other) const;
    bool operator!=(const BennettsFreeEnergyAverage &other) const;
    
    BennettsFreeEnergyAverage operator+(const BennettsFreeEnergyAverage &other) const;
    
    BennettsFreeEnergyAverage& operator+=(const BennettsFreeEnergyAverage &other);
    
    static const char* typeName();

    QString toString() const;

    void clear();
    
    void accumulate(double value);

    double forwardsRatio() const;
    double forwardsStandardError(double level) const;
    
    double backwardsRatio() const;
    double backwardsStandardError(double level) const;

private:
    /** The average forwards ratio */
    double fwds_avg;
    
    /** The average backwards ratio */
    double bwds_avg;
    
    /** The average of the squared forwards ratio */
    double fwds_avg2;
    
    /** The average of the squared backwards ratio */
    double bwds_avg2;
};

}

Q_DECLARE_METATYPE( SireMaths::FreeEnergyAverage )
Q_DECLARE_METATYPE( SireMaths::BennettsFreeEnergyAverage )

SIRE_EXPOSE_CLASS( SireMaths::FreeEnergyAverage )
SIRE_EXPOSE_CLASS( SireMaths::BennettsFreeEnergyAverage )

SIRE_END_HEADER

#endif
