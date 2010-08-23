#ifndef SIRENGLOBAL_H
#define SIRENGLOBAL_H

/** This file contains some global macros that are used
    to improve the portability of the code. */

//include qglobal.h - this includes all of the definitions
//that are needed for these macros
#include <qconfig.h>
#include <qglobal.h>

#ifdef __cplusplus

#include <qmetatype.h>

#include <boost/current_function.hpp>
//macro used to return the current file and line as a QString
#ifdef CODELOC
#undef CODELOC
#endif

#define CODELOC \
   QObject::tr("FILE: %1, LINE: %2, FUNCTION: %3") \
           .arg(__FILE__).arg(__LINE__).arg(BOOST_CURRENT_FUNCTION)

#ifdef QUICK_CODELOC
#undef QUICK_CODELOC
#endif

#define QUICK_CODELOC __FILE__, __LINE__, BOOST_CURRENT_FUNCTION

#ifdef QUICK_CODELOC_ARGS
#undef QUICK_CODELOC_ARGS
#endif

#define QUICK_CODELOC_ARGS const char *file, int line, const char *current_function

//copy the QT_BEGIN_HEADER and QT_END_HEADER
//to SIRE_BEGIN_HEADER and SIRE_END_HEADER. This
//will allow me to change their definition at some
//future point in time without having to mess with Qt
#define SIREN_BEGIN_HEADER QT_BEGIN_HEADER
#define SIREN_END_HEADER QT_END_HEADER

//create keywords that control whether classes or functions are
//exposed to a scripting language
#define SIREN_EXPOSE_FUNCTION(f)  /* Exposing function #1 */
#define SIREN_EXPOSE_CLASS(c)     /* Exposing class #1 */
#define SIREN_EXPOSE_ALIAS(c,a)   /* Exposing class #1 as alias #2 */
#define SIREN_EXPOSE_OBJECT_PTR(c,a)  /* Exposing pointer #1 to object class #2 */
#define SIREN_EXPOSE_HANDLE_PTR(c,a)  /* Exposing pointer #1 to object class #2 */

#ifdef SIREN_INSTANTIATE_TEMPLATES
    #define SIREN_EXPOSE_PRIMITIVE(X) \
        namespace Siren { \
        template class Siren::Primitive<X>; template class Siren::PrimitiveObject<X>; } \
        Q_DECLARE_METATYPE( Siren::PrimitiveObject< X > )

    #define SIREN_EXPOSE_ALIAS_PRIMITIVE(X,Y) \
        namespace Siren { \
        template class Siren::Primitive<Y>; template class Siren::PrimitiveObject<Y>; } \
        Q_DECLARE_METATYPE( Siren::PrimitiveObject<Y> )
#else
    #define SIREN_EXPOSE_PRIMITIVE(X) Q_DECLARE_METATYPE( Siren::PrimitiveObject< X > )
    #define SIREN_EXPOSE_ALIAS_PRIMITIVE(X,Y) \
        Q_DECLARE_METATYPE( Siren::PrimitiveObject<Y> )
#endif

//create the keyword used to export a symbol - this
//is a copy of Q_DECL_EXPORT, which will definitely
//be set to the correct value
#ifndef SIREN_NO_VISIBILITY_AVAILABLE
#define SIREN_EXPORT Q_DECL_EXPORT
#define SIRENTEST_EXPORT Q_DECL_EXPORT

//create the keyword to import a symbol - this is copied
//from Q_DECL_IMPORT, which will definitely be set to
//the correct value
#define SIREN_IMPORT Q_DECL_IMPORT
#define SIRENTEST_IMPORT Q_DECL_IMPORT

#else
#define SIREN_EXPORT
#define SIRENTEST_EXPORT
#define SIREN_IMPORT
#define SIRENTEST_IMPORT
#endif

//create the keyword to fix symbol visibility problems for out-of-line
//template functions
#define SIREN_OUTOFLINE_TEMPLATE Q_OUTOFLINE_TEMPLATE

//do the same for inline template functions
#define SIREN_INLINE_TEMPLATE Q_INLINE_TEMPLATE

//Sire is much more picky about what the compiler can support
//than Qt - use the results of the Qt portability tests
//to check that the C++ compiler is up to scratch

#ifdef Q_BROKEN_TEMPLATE_SPECIALIZATION
#error Siren requires a C++ compiler that can handle templates in \
       a standards compliant manner. Please use a more recent version \
       of your compiler, or try gcc >= 3.0
#endif

#ifdef QT_NO_PARTIAL_TEMPLATE_SPECIALIZATION
#error Siren requires that the C++ compiler supports \
       partial template specialization. Please use a more \
     recent version of your compiler, or try gcc >= 3.0
#endif

#ifdef Q_NO_USING_KEYWORD
#error Siren requires that the C++ compiler supports the 'using' \
       keyword. Please use a more recent version of your compiler \
     or try gcc >= 3.0
#endif

#ifdef QT_NO_MEMBER_TEMPLATES
#error Siren requires that the C++ compiler supports the use of \
       member templates. Please use a more recent version of your \
     compiler or try gcc >= 3.0
#endif

#ifdef QT_NO_THREAD
#error Siren requires that Qt was built with threading enabled. \
       Please recompile Qt with threading enabled.
#endif

#ifdef QT_NO_TEXTSTREAM
#error Siren requires that Qt was built with QTextStream enabled \
       Please recompile Qt with QTextStream enabled.
#endif

//sire requires that at least the QtCore and QtSql modules
//are available
#if  !(QT_EDITION & QT_MODULE_CORE)
#error Sire requires that the QtCore module is enabled.
#endif

namespace Siren
{

/** These functions convert a pointer into an integer.

    toInt(void*) : Returns an integer from the pointer. Returns
                   an integer of the same size as the pointer.

    toUInt(void*) : Returns an unsigned integer of the same size as
                    the pointer.

    toInt32(void*)   : Returns a 32bit integer from a pointer. This
                       will truncate if the pointer is larger than 32bits.

    toUInt32(void*)  : Same as toInt32, but returns the pointer as an
                       unsigned 32bit integer.

    toInt64(void*) : Returns a 64bit integer

    toUInt64(void*) : Returns a 64bit unsigned integer.

    @author Christopher Woods
*/
#ifdef QT_POINTER_SIZE
  #if QT_POINTER_SIZE == 4

    // Functions used if we have 32bit pointers
    inline qint32 toInt(const void *ptr)
    {
        return qint32(ptr);
    }

    inline quint32 toUInt(const void *ptr)
    {
        return quint32(ptr);
    }

    inline qint32 toInt32(const void *ptr)
    {
        return toInt(ptr);
    }

    inline quint32 toUInt32(const void *ptr)
    {
        return toUInt(ptr);
    }

    inline qint64 toInt64(const void *ptr)
    {
        return qint64( toInt(ptr) );
    }

    inline quint64 toUInt64(const void *ptr)
    {
        return quint64( toUInt(ptr) );
    }

  #elif QT_POINTER_SIZE == 8

    // Functions used if we have 64bit pointers
    inline qint64 toInt(const void *ptr)
    {
        return qint64(ptr);
    }

    inline quint64 toUInt(const void *ptr)
    {
        return quint64(ptr);
    }

    inline qint32 toInt32(const void *ptr)
    {
        return qint32( toInt(ptr) );
    }

    inline quint32 toUInt32(const void *ptr)
    {
        return qint32( toUInt(ptr) );
    }

    inline qint64 toInt64(const void *ptr)
    {
        return toInt(ptr);
    }

    inline quint64 toUInt64(const void *ptr)
    {
        return toUInt(ptr);
    }

  #else
    #fatal Invalid QT_POINTER_SIZE value stored!
  #endif
#else
  #fatal No QT_POINTER_SIZE macro defined!
#endif

} // end of namespace Siren

//Add functions that are used to register all public
//types (this allows them to be streamed to a binary file and
//created dynamically)
#include <QMetaType>

namespace Siren
{

/** This is the base class of all metatype registrations */
class RegisterMetaType
{
protected:
    enum { IS_OBJECT = 1, IS_HANDLE = 2, IS_PRIMITIVE = 3 };

public:
    ~RegisterMetaType()
    {}

    const QString& typeName() const
    {
        return type_name;
    }

    bool isObject() const
    {
        return sub_type == IS_OBJECT;
    }

    bool isHandle() const
    {
        return sub_type == IS_HANDLE;
    }

    bool isPrimitive() const
    {
        return sub_type == IS_PRIMITIVE;
    }

    bool isConcrete() const
    {
        return is_concrete;
    }

    virtual RegisterMetaType* clone() const=0;

protected:
    RegisterMetaType(int subtype, const QString &name, 
                     bool reg=true) : type_name(name), sub_type(subtype)
    {
        if (reg)
            is_concrete = true;
        else
            is_concrete = false;
    }
    
    RegisterMetaType(const RegisterMetaType &other)
          : type_name(other.type_name), sub_type(other.sub_type),
            is_concrete(other.is_concrete)
    {}

    static const RegisterMetaType* getRegistration(const QString &classname);

    static void registerClassName(const QString &name, 
                                  const RegisterMetaType *metatype);

private:
    QString type_name;
    qint8 sub_type;
    bool is_concrete;
};


namespace detail
{
struct VIRTUAL_CLASS_TYPE{ };
struct NONSTREAMING_CLASS_TYPE{ };
} // end of namespace detail

static const detail::VIRTUAL_CLASS_TYPE VIRTUAL_CLASS = detail::VIRTUAL_CLASS_TYPE();
static const detail::NONSTREAMING_CLASS_TYPE NONSTREAMING_CLASS = 
                                                detail::NONSTREAMING_CLASS_TYPE();

/** This is used to register a Siren::Object of type 'T' - this
    needs to be called once for each public Siren::Object

    @author Christopher Woods
*/
template<class T>
class RegisterObject : public RegisterMetaType
{
public:
    /** Use this constructor to register a concrete Object */
    RegisterObject()
        : RegisterMetaType( 
                    RegisterMetaType::IS_OBJECT,
                    QMetaType::typeName( qMetaTypeId<T>() ) )
    {
        qRegisterMetaType<T>(this->typeName().toAscii().constData());
        qRegisterMetaTypeStreamOperators<T>(this->typeName().toAscii().constData());

        RegisterMetaType::registerClassName(this->typeName(), this);
    }

    /** Use this constructor to register a virtual Object */
    RegisterObject( const detail::VIRTUAL_CLASS_TYPE& )
        : RegisterMetaType( RegisterMetaType::IS_OBJECT, T::typeName(), false )
    {
        RegisterMetaType::registerClassName(this->typeName(), this);
    }

    RegisterObject(const RegisterObject<T> &other)
        : RegisterMetaType(other)
    {}

    /** Destructor */
    ~RegisterObject()
    {}

    RegisterObject<T>* clone() const
    {
        return new RegisterObject<T>(*this);
    }

    static const RegisterMetaType* getRegistration()
    {
        return RegisterMetaType::getRegistration( T::typeName() );
    }
};

/** This is used to register a Siren::Handle of type 'T' - this
    needs to be called once for each public Siren::Handle

    @author Christopher Woods
*/
template<class T>
class RegisterHandle : public RegisterMetaType
{
public:
    /** Use this constructor to register a concrete Handle */
    RegisterHandle()
        : RegisterMetaType( 
                    RegisterMetaType::IS_HANDLE,
                    QMetaType::typeName( qMetaTypeId<T>() ) )
    {
        qRegisterMetaType<T>(this->typeName().toAscii().constData());
        RegisterMetaType::registerClassName(this->typeName(), this);
    }

    /** Use this constructor to register a virtual Handle */
    RegisterHandle( const detail::VIRTUAL_CLASS_TYPE& )
        : RegisterMetaType( RegisterMetaType::IS_HANDLE, T::typeName(), false )
    {
        RegisterMetaType::registerClassName(this->typeName(), this);
    }
    
    RegisterHandle(const RegisterHandle<T> &other) : RegisterMetaType(other)
    {}

    /** Destructor */
    ~RegisterHandle()
    {}

    RegisterHandle<T>* clone() const
    {
        return new RegisterHandle<T>(*this);
    }

    static const RegisterMetaType* getRegistration()
    {
        return RegisterMetaType::getRegistration( T::typeName() );
    }
};

template<class T> class PrimitiveObject;

/** This is the base class of RegisterPrimitive<T> */
class RegisterPrimitiveBase : public RegisterMetaType
{
public:
    ~RegisterPrimitiveBase()
    {}
    
    const RegisterMetaType* getObjectRegistration() const
    {
        return primitive_object;
    }
    
protected:
    RegisterPrimitiveBase(int subtype, const QString &name)
            : RegisterMetaType(subtype, name, true)
    {}

    RegisterPrimitiveBase(const RegisterPrimitiveBase &other)
            : RegisterMetaType(other), primitive_object(other.primitive_object)
    {}

    const RegisterMetaType *primitive_object;
};

/** This is used to register a Siren::Primitive<T> derived type 
    of type 'T' - this needs to be called once for each public Siren::Primitive

    @author Christopher Woods
*/
template<class T>
class RegisterPrimitive : public RegisterPrimitiveBase
{
public:
    /** Use this constructor to register a Primitive */
    RegisterPrimitive()
        : RegisterPrimitiveBase( 
                    RegisterMetaType::IS_PRIMITIVE,
                    QMetaType::typeName( qMetaTypeId<T>() ) )
    {
        qRegisterMetaType<T>(this->typeName().toAscii().constData());
        qRegisterMetaTypeStreamOperators<T>(this->typeName().toAscii().constData());
        
        RegisterObject< PrimitiveObject<T> > reg_object;
        primitive_object = RegisterObject< PrimitiveObject<T> >::getRegistration();

        RegisterMetaType::registerClassName(this->typeName(), this);
    }

    RegisterPrimitive(const RegisterPrimitive<T> &other)
         : RegisterPrimitiveBase(other)
    {}

    /** Destructor */
    ~RegisterPrimitive()
    {}

    RegisterPrimitive<T>* clone() const
    {
        return new RegisterPrimitive<T>(*this);
    }

    static const RegisterMetaType* getRegistration()
    {
        return RegisterMetaType::getRegistration( QMetaType::typeName(qMetaTypeId<T>()) );
    }
};

} // end of namespace Siren

#else  // else #ifdef __cplusplus

//copied directly from qglobal.h /////////////////
#ifndef Q_DECL_EXPORT
#  ifdef Q_OS_WIN
#    define Q_DECL_EXPORT __declspec(dllexport)
#  elif defined(QT_VISIBILITY_AVAILABLE)
#    define Q_DECL_EXPORT __attribute__((visibility("default")))
#  endif
#  ifndef Q_DECL_EXPORT
#    define Q_DECL_EXPORT
#  endif
#endif
#ifndef Q_DECL_IMPORT
#  ifdef Q_OS_WIN
#    define Q_DECL_IMPORT __declspec(dllimport)
#  else
#    define Q_DECL_IMPORT
#  endif
#endif
///////////////////////////////////////////////////

#ifndef SIREN_NO_VISIBILITY_AVAILABLE
#define SIREN_EXPORT Q_DECL_EXPORT
#define SIREN_IMPORT Q_DECL_IMPORT
#define SIRENTEST_EXPORT Q_DECL_EXPORT
#define SIRENTEST_IMPORT Q_DECL_IMPORT
#else
#define SIREN_EXPORT
#define SIREN_IMPORT
#define SIRENTEST_EXPORT
#define SIRENTEST_IMPORT
#endif

#endif // #ifdef __cplusplus

#endif // SIREGLOBAL_H
