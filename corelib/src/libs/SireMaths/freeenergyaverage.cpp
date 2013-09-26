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

#include "freeenergyaverage.h"

#include "SireUnits/units.h"

#include "SireMaths/maths.h"
#include "SireMaths/histogram.h"

#include "SireStream/datastream.h"

#include <QDebug>

using namespace SireMaths;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireBase;
using namespace SireStream;

////////////
//////////// Implementation of FreeEnergyAverage
////////////

static const RegisterMetaType<FreeEnergyAverage> r_avg;

/** Serialise to a binary datastream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds,
                                         const FreeEnergyAverage &avg)
{
    writeHeader(ds, r_avg, 2);
    
    SharedDataStream sds(ds);
    
    sds << avg.hist << static_cast<const ExpAverage&>(avg);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, FreeEnergyAverage &avg)
{
    VersionID v = readHeader(ds, r_avg);
    
    if (v == 2)
    {
        SharedDataStream sds(ds);
        
        sds >> avg.hist >> static_cast<ExpAverage&>(avg);
    }
    else if (v == 1)
    {
        ds >> static_cast<ExpAverage&>(avg);
        avg.hist = Histogram();
    }
    else
        throw version_error(v, "1", r_avg, CODELOC);
        
    return ds;
}

/** Constructor - this defaults to accumulating the average
    at room temperature (25 C) and collects statistics about the
    free energy using a histogram of bin width 0.5 kcal mol-1 */
FreeEnergyAverage::FreeEnergyAverage()
                  : ConcreteProperty<FreeEnergyAverage,ExpAverage>(
                        -1.0 / (k_boltz * double(25*celsius)) ), hist(0.5)
{}

/** Construct an accumulator to accumulate the free energy average
    at the specified temperature, and to collect statistics about
    the free energy using a histogram of bin width 0.5 kcal mol-1 */
FreeEnergyAverage::FreeEnergyAverage(const Temperature &temperature)
                  : ConcreteProperty<FreeEnergyAverage,ExpAverage>(
                        -1.0 / (k_boltz * temperature.to(kelvin)) ), hist(0.5)
{}

/** Constructor - this defaults to accumulating the average
    at room temperature (25 C) and collects statistics about the
    free energy using a histogram of the passed bin width. If the binwidth
    is zero, then a histogram of energies is not collected */
FreeEnergyAverage::FreeEnergyAverage(const MolarEnergy &binwidth)
                  : ConcreteProperty<FreeEnergyAverage,ExpAverage>(
                        -1.0 / (k_boltz * double(25*celsius)) ), hist(binwidth.value())
{}

/** Construct an accumulator to accumulate the free energy average
    at the specified temperature, and to collect statistics about
    the free energy using a histogram of passed bin width */
FreeEnergyAverage::FreeEnergyAverage(const Temperature &temperature,
                                     const MolarEnergy &binwidth)
                  : ConcreteProperty<FreeEnergyAverage,ExpAverage>(
                        -1.0 / (k_boltz * temperature.to(kelvin)) ), hist(binwidth.value())
{}

/** Copy constructor */
FreeEnergyAverage::FreeEnergyAverage(const FreeEnergyAverage &other)
                  : ConcreteProperty<FreeEnergyAverage,ExpAverage>(other),
                    hist(other.hist)
{}

/** Destructor */
FreeEnergyAverage::~FreeEnergyAverage()
{}

/** Copy assignment operator */
FreeEnergyAverage& FreeEnergyAverage::operator=(const FreeEnergyAverage &other)
{
    if (this != &other)
    {
        hist = other.hist;
        ExpAverage::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool FreeEnergyAverage::operator==(const FreeEnergyAverage &other) const
{
    return this == &other or
           (hist == other.hist and ExpAverage::operator==(other));
}

/** Comparison operator */
bool FreeEnergyAverage::operator!=(const FreeEnergyAverage &other) const
{
    return not operator==(other);
}

/** Combine the passed average onto this average */
FreeEnergyAverage& FreeEnergyAverage::operator+=(const FreeEnergyAverage &other)
{
    ExpAverage::operator+=(other);
    hist += other.hist;
    return *this;
}

/** Return the combination of this average plus other */
FreeEnergyAverage FreeEnergyAverage::operator+(const FreeEnergyAverage &other) const
{
    FreeEnergyAverage ret(*this);
    ret += other;
    return ret;
}

/** Return the temperature at which the free energy average
    is being accumulated */
Temperature FreeEnergyAverage::temperature() const
{
    return Temperature( -1.0 / (k_boltz*scaleFactor()) );
}

/** Return the histogram of energies */
const Histogram& FreeEnergyAverage::histogram() const
{
    return hist;
}

QString FreeEnergyAverage::toString() const
{
    return QObject::tr("FreeEnergyAverage( dG = %1 kcal mol-1, average = %2 kcal mol-1 "
                       "stderr = %3 kcal mol-1, "
                       "skew = %4 kcal mol-1, nSamples = %5 )")
                            .arg(this->average())
                            .arg(histogram().mean())
                            .arg(histogram().standardDeviation())
                            .arg(histogram().skew())
                            .arg(nSamples());
}

const char* FreeEnergyAverage::typeName()
{
    return QMetaType::typeName( qMetaTypeId<FreeEnergyAverage>() );
}

/** Clear all data from the accumulator */
void FreeEnergyAverage::clear()
{
    hist = Histogram( hist.binWidth() );
    ExpAverage::clear();
}

/** Accumulate the passed energy difference onto the free energy average */
void FreeEnergyAverage::accumulate(double value)
{
    ExpAverage::accumulate(value);
    hist.accumulate(value);
}

/** Return the Taylor series expansion estimate the difference in free energy */
double FreeEnergyAverage::taylorExpansion() const
{
    return hist.mean() - 0.5*k_boltz*temperature() *
              ( hist.meanOfSquares() - (hist.mean()*hist.mean()) );
}

////////////
//////////// Implementation of BennettsFreeEnergyAverage
////////////

static const RegisterMetaType<BennettsFreeEnergyAverage> r_bennetts;

QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds,
                                         const BennettsFreeEnergyAverage &bennetts)
{
    writeHeader(ds, r_bennetts, 1);
    
    SharedDataStream sds(ds);
    
    sds << bennetts.fwds_avg << bennetts.bwds_avg
        << bennetts.fwds_avg2 << bennetts.bwds_avg2
        << static_cast<const FreeEnergyAverage&>(bennetts);
    
    return ds;
}

QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, BennettsFreeEnergyAverage &bennetts)
{
    VersionID v = readHeader(ds, r_bennetts);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> bennetts.fwds_avg >> bennetts.bwds_avg
            >> bennetts.fwds_avg2 >> bennetts.bwds_avg2
            >> static_cast<FreeEnergyAverage&>(bennetts);
    }
    else
        throw version_error(v, "1", r_bennetts, CODELOC);
    
    return ds;
}

/** Constructor */
BennettsFreeEnergyAverage::BennettsFreeEnergyAverage()
     : ConcreteProperty<BennettsFreeEnergyAverage,FreeEnergyAverage>(),
       fwds_avg(0), bwds_avg(0), fwds_avg2(0), bwds_avg2(0)
{}

/** Construct the average at the specified temperature */
BennettsFreeEnergyAverage::BennettsFreeEnergyAverage(const Temperature &temperature)
     : ConcreteProperty<BennettsFreeEnergyAverage,FreeEnergyAverage>(temperature),
       fwds_avg(0), bwds_avg(0), fwds_avg2(0), bwds_avg2(0)
{}

/** Construct the average using a histogram of the specified bin width */
BennettsFreeEnergyAverage::BennettsFreeEnergyAverage(const MolarEnergy &binwidth)
     : ConcreteProperty<BennettsFreeEnergyAverage,FreeEnergyAverage>(binwidth),
       fwds_avg(0), bwds_avg(0), fwds_avg2(0), bwds_avg2(0)
{}

/** Construct at the specified temperature, using a histogram of the specified bin width */
BennettsFreeEnergyAverage::BennettsFreeEnergyAverage(const Temperature &temperature,
                                                   const MolarEnergy &binwidth)
     : ConcreteProperty<BennettsFreeEnergyAverage,FreeEnergyAverage>(temperature,binwidth),
       fwds_avg(0), bwds_avg(0), fwds_avg2(0), bwds_avg2(0)
{}

/** Copy constructor */
BennettsFreeEnergyAverage::BennettsFreeEnergyAverage(const BennettsFreeEnergyAverage &other)
     : ConcreteProperty<BennettsFreeEnergyAverage,FreeEnergyAverage>(other),
       fwds_avg(other.fwds_avg), bwds_avg(other.bwds_avg),
       fwds_avg2(other.fwds_avg2), bwds_avg2(other.bwds_avg2)
{}

/** Destructor */
BennettsFreeEnergyAverage::~BennettsFreeEnergyAverage()
{}

/** Copy assignment operator */
BennettsFreeEnergyAverage&
BennettsFreeEnergyAverage::operator=(const BennettsFreeEnergyAverage &other)
{
    if (this != &other)
    {
        fwds_avg = other.fwds_avg;
        bwds_avg = other.bwds_avg;
        fwds_avg2 = other.fwds_avg2;
        bwds_avg2 = other.bwds_avg2;
        FreeEnergyAverage::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool BennettsFreeEnergyAverage::operator==(const BennettsFreeEnergyAverage &other) const
{
    return fwds_avg == other.fwds_avg and bwds_avg == other.bwds_avg and
           fwds_avg2 == other.fwds_avg2 and bwds_avg2 == other.bwds_avg2 and
           FreeEnergyAverage::operator==(other);
}

/** Comparison operator */
bool BennettsFreeEnergyAverage::operator!=(const BennettsFreeEnergyAverage &other) const
{
    return not operator==(other);
}

/** Self-addition operator */
BennettsFreeEnergyAverage&
BennettsFreeEnergyAverage::operator+=(const BennettsFreeEnergyAverage &other)
{
    double nsteps = nSamples() + other.nSamples();
        
    double my_ratio = nSamples() / nsteps;
    double other_ratio = other.nSamples() / nsteps;

    FreeEnergyAverage::operator+=(other);
    
    fwds_avg = fwds_avg * my_ratio + other.fwds_avg * other_ratio;
    bwds_avg = bwds_avg * my_ratio + other.bwds_avg * other_ratio;
    fwds_avg2 = fwds_avg2 * my_ratio + other.fwds_avg2 * other_ratio;
    bwds_avg2 = bwds_avg2 * my_ratio + other.bwds_avg2 * other_ratio;
    
    return *this;
}

/** Addition operator */
BennettsFreeEnergyAverage
BennettsFreeEnergyAverage::operator+(const BennettsFreeEnergyAverage &other) const
{
    BennettsFreeEnergyAverage ret(*this);
    ret += other;
    return ret;
}

const char* BennettsFreeEnergyAverage::typeName()
{
    return QMetaType::typeName( qMetaTypeId<BennettsFreeEnergyAverage>() );
}

QString BennettsFreeEnergyAverage::toString() const
{
    return QObject::tr("BennettsFreeEnergyAverage( dG = %1 kcal mol-1, average = %2 kcal mol-1 "
                       "forwardsRatio() = %6, backwardsRatio() = %7, stderr = %3 kcal mol-1, "
                       "skew = %4 kcal mol-1, nSamples = %5 )")
                            .arg(this->average())
                            .arg(histogram().mean())
                            .arg(histogram().standardDeviation())
                            .arg(histogram().skew())
                            .arg(nSamples())
                            .arg(forwardsRatio())
                            .arg(backwardsRatio());
}

/** Clear this accumulator */
void BennettsFreeEnergyAverage::clear()
{
    FreeEnergyAverage::clear();
    fwds_avg = 0;
    bwds_avg = 0;
    fwds_avg2 = 0;
    bwds_avg2 = 0;
}

/** Accumulate the passed value onto the average */
void BennettsFreeEnergyAverage::accumulate(double value)
{
    double nsteps = nSamples() + 1;
        
    double my_ratio = nSamples() / nsteps;
    double other_ratio = 1.0 / nsteps;

    double fwds_val = 1.0 / (1.0 + std::exp(this->scaleFactor()*value));
    double bwds_val = 1.0 / (1.0 + std::exp(-this->scaleFactor()*value));
    
    fwds_avg = my_ratio*fwds_avg + other_ratio*fwds_val;
    bwds_avg = my_ratio*bwds_avg + other_ratio*bwds_val;
    
    fwds_avg2 = my_ratio*fwds_avg2 + other_ratio*(fwds_val*fwds_val);
    bwds_avg2 = my_ratio*bwds_avg2 + other_ratio*(bwds_val*bwds_val);
    
    FreeEnergyAverage::accumulate(value);
}

/** Return the forwards part of the ratio. This is the ensemble average
    of 1 / (1 + exp(-dE/kT)) */
double BennettsFreeEnergyAverage::forwardsRatio() const
{
    return fwds_avg;
}

/** Return the standard error on the forwards value to the passed confidence level */
double BennettsFreeEnergyAverage::forwardsStandardError(double level) const
{
    if (this->nSamples() == 0)
        return 0;

    double stdev = std::sqrt( fwds_avg2 - pow_2(fwds_avg) );
    return Histogram::tValue(this->nSamples(),level) * stdev / std::sqrt(this->nSamples());
}

/** Return the backwards part of the ratio. This is the ensemble average
    of 1 / (1 + exp(dE/kT)) */
double BennettsFreeEnergyAverage::backwardsRatio() const
{
    return bwds_avg;
}

/** Return the standard error on the backwards value to the passed confidence level */
double BennettsFreeEnergyAverage::backwardsStandardError(double level) const
{
    if (this->nSamples() == 0)
        return 0;

    double stdev = std::sqrt( bwds_avg2 - pow_2(bwds_avg) );
    return Histogram::tValue(this->nSamples(),level) * stdev / std::sqrt(this->nSamples());
}
