
#include "i.h"
#include "registerexpression.h"
#include "expression.h"
#include "complexvalues.h"

#include "SireMaths/errors.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireCAS;

static const RegisterMetaType<SireCAS::I> r_i;

/** Serialise to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const I&)
{
    writeHeader(ds, r_i, 0);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, I&)
{
    VersionID v = readHeader(ds, r_i);

    if (v != 0)
        throw version_error(v, "0", r_i, CODELOC);

    return ds;
}

/** Register this expression */
static RegisterExpression<I> RegisterI;

/** Constructor */
I::I() : Constant()
{}

/** Copy constructor */
I::I(const I &other) : Constant(other)
{}

/** Destructor */
I::~I()
{}

/** Comparison operator */
bool I::operator==(const ExBase &other) const
{
    const I *other_i = dynamic_cast<const I*>(&other);

    return other_i != 0 and typeid(other).name() == typeid(*this).name();
}

/** Return a hash of this expression */
uint I::hash() const
{
    return ( r_i.magicID() << 16 ) | ( r_i.magicID() & 0x0000FFFF );
}

/** Return a string representation */
QString I::toString() const
{
    return "i";
}

/** Cannot evaluate 'i' as a real number, so throw a domain error */
double I::evaluate(const Values&) const
{
    throw SireMaths::domain_error(QObject::tr(
        "Cannot evaluate 'i' as a real number"), CODELOC);
}

/** Evaluate this as a complex number - return 'i' */
Complex I::evaluate(const ComplexValues&) const
{
    return Complex(0,1);
}

/** Return the complex conjugate of 'i' (-i) */
Expression I::conjugate() const
{
    return -1 * toExpression();
}

/** I is definitely complex :-) */
bool I::isComplex() const
{
    return true;
}
