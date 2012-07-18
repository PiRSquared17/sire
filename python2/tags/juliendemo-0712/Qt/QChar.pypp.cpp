// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include <QString>
#include <QByteArray>
#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QTextStream>
#include <QDateTime>
#include <QLocale>
#include <QUuid>
#include <qnamespace.h>
#include <QVariant>
#include <QUrl>
#include <QBitArray>
#include "QChar.pypp.hpp"

namespace bp = boost::python;

PyObject* convert(QChar const& c);

QChar __copy__(const QChar &other){ return QChar(other); }

#include "Qt/qdatastream.hpp"

void register_QChar_class(){

    { //::QChar
        typedef bp::class_< QChar > QChar_exposer_t;
        QChar_exposer_t QChar_exposer = QChar_exposer_t( "QChar" );
        bp::scope QChar_scope( QChar_exposer );
        bp::enum_< QChar::Category>("Category")
            .value("NoCategory", QChar::NoCategory)
            .value("Mark_NonSpacing", QChar::Mark_NonSpacing)
            .value("Mark_SpacingCombining", QChar::Mark_SpacingCombining)
            .value("Mark_Enclosing", QChar::Mark_Enclosing)
            .value("Number_DecimalDigit", QChar::Number_DecimalDigit)
            .value("Number_Letter", QChar::Number_Letter)
            .value("Number_Other", QChar::Number_Other)
            .value("Separator_Space", QChar::Separator_Space)
            .value("Separator_Line", QChar::Separator_Line)
            .value("Separator_Paragraph", QChar::Separator_Paragraph)
            .value("Other_Control", QChar::Other_Control)
            .value("Other_Format", QChar::Other_Format)
            .value("Other_Surrogate", QChar::Other_Surrogate)
            .value("Other_PrivateUse", QChar::Other_PrivateUse)
            .value("Other_NotAssigned", QChar::Other_NotAssigned)
            .value("Letter_Uppercase", QChar::Letter_Uppercase)
            .value("Letter_Lowercase", QChar::Letter_Lowercase)
            .value("Letter_Titlecase", QChar::Letter_Titlecase)
            .value("Letter_Modifier", QChar::Letter_Modifier)
            .value("Letter_Other", QChar::Letter_Other)
            .value("Punctuation_Connector", QChar::Punctuation_Connector)
            .value("Punctuation_Dash", QChar::Punctuation_Dash)
            .value("Punctuation_Open", QChar::Punctuation_Open)
            .value("Punctuation_Close", QChar::Punctuation_Close)
            .value("Punctuation_InitialQuote", QChar::Punctuation_InitialQuote)
            .value("Punctuation_FinalQuote", QChar::Punctuation_FinalQuote)
            .value("Punctuation_Other", QChar::Punctuation_Other)
            .value("Symbol_Math", QChar::Symbol_Math)
            .value("Symbol_Currency", QChar::Symbol_Currency)
            .value("Symbol_Modifier", QChar::Symbol_Modifier)
            .value("Symbol_Other", QChar::Symbol_Other)
            .value("Punctuation_Dask", QChar::Punctuation_Dask)
            .export_values()
            ;
        bp::enum_< QChar::CombiningClass>("CombiningClass")
            .value("Combining_BelowLeftAttached", QChar::Combining_BelowLeftAttached)
            .value("Combining_BelowAttached", QChar::Combining_BelowAttached)
            .value("Combining_BelowRightAttached", QChar::Combining_BelowRightAttached)
            .value("Combining_LeftAttached", QChar::Combining_LeftAttached)
            .value("Combining_RightAttached", QChar::Combining_RightAttached)
            .value("Combining_AboveLeftAttached", QChar::Combining_AboveLeftAttached)
            .value("Combining_AboveAttached", QChar::Combining_AboveAttached)
            .value("Combining_AboveRightAttached", QChar::Combining_AboveRightAttached)
            .value("Combining_BelowLeft", QChar::Combining_BelowLeft)
            .value("Combining_Below", QChar::Combining_Below)
            .value("Combining_BelowRight", QChar::Combining_BelowRight)
            .value("Combining_Left", QChar::Combining_Left)
            .value("Combining_Right", QChar::Combining_Right)
            .value("Combining_AboveLeft", QChar::Combining_AboveLeft)
            .value("Combining_Above", QChar::Combining_Above)
            .value("Combining_AboveRight", QChar::Combining_AboveRight)
            .value("Combining_DoubleBelow", QChar::Combining_DoubleBelow)
            .value("Combining_DoubleAbove", QChar::Combining_DoubleAbove)
            .value("Combining_IotaSubscript", QChar::Combining_IotaSubscript)
            .export_values()
            ;
        bp::enum_< QChar::Decomposition>("Decomposition")
            .value("NoDecomposition", QChar::NoDecomposition)
            .value("Canonical", QChar::Canonical)
            .value("Font", QChar::Font)
            .value("NoBreak", QChar::NoBreak)
            .value("Initial", QChar::Initial)
            .value("Medial", QChar::Medial)
            .value("Final", QChar::Final)
            .value("Isolated", QChar::Isolated)
            .value("Circle", QChar::Circle)
            .value("Super", QChar::Super)
            .value("Sub", QChar::Sub)
            .value("Vertical", QChar::Vertical)
            .value("Wide", QChar::Wide)
            .value("Narrow", QChar::Narrow)
            .value("Small", QChar::Small)
            .value("Square", QChar::Square)
            .value("Compat", QChar::Compat)
            .value("Fraction", QChar::Fraction)
            .export_values()
            ;
        bp::enum_< QChar::Direction>("Direction")
            .value("DirL", QChar::DirL)
            .value("DirR", QChar::DirR)
            .value("DirEN", QChar::DirEN)
            .value("DirES", QChar::DirES)
            .value("DirET", QChar::DirET)
            .value("DirAN", QChar::DirAN)
            .value("DirCS", QChar::DirCS)
            .value("DirB", QChar::DirB)
            .value("DirS", QChar::DirS)
            .value("DirWS", QChar::DirWS)
            .value("DirON", QChar::DirON)
            .value("DirLRE", QChar::DirLRE)
            .value("DirLRO", QChar::DirLRO)
            .value("DirAL", QChar::DirAL)
            .value("DirRLE", QChar::DirRLE)
            .value("DirRLO", QChar::DirRLO)
            .value("DirPDF", QChar::DirPDF)
            .value("DirNSM", QChar::DirNSM)
            .value("DirBN", QChar::DirBN)
            .export_values()
            ;
        bp::enum_< QChar::Joining>("Joining")
            .value("OtherJoining", QChar::OtherJoining)
            .value("Dual", QChar::Dual)
            .value("Right", QChar::Right)
            .value("Center", QChar::Center)
            .export_values()
            ;
        bp::enum_< QChar::SpecialCharacter>("SpecialCharacter")
            .value("Null", QChar::Null)
            .value("Nbsp", QChar::Nbsp)
            .value("ReplacementCharacter", QChar::ReplacementCharacter)
            .value("ObjectReplacementCharacter", QChar::ObjectReplacementCharacter)
            .value("ByteOrderMark", QChar::ByteOrderMark)
            .value("ByteOrderSwapped", QChar::ByteOrderSwapped)
            .value("ParagraphSeparator", QChar::ParagraphSeparator)
            .value("LineSeparator", QChar::LineSeparator)
            .export_values()
            ;
        bp::enum_< QChar::UnicodeVersion>("UnicodeVersion")
            .value("Unicode_Unassigned", QChar::Unicode_Unassigned)
            .value("Unicode_1_1", QChar::Unicode_1_1)
            .value("Unicode_2_0", QChar::Unicode_2_0)
            .value("Unicode_2_1_2", QChar::Unicode_2_1_2)
            .value("Unicode_3_0", QChar::Unicode_3_0)
            .value("Unicode_3_1", QChar::Unicode_3_1)
            .value("Unicode_3_2", QChar::Unicode_3_2)
            .value("Unicode_4_0", QChar::Unicode_4_0)
            .export_values()
            ;
        QChar_exposer.def( bp::init< >() );
        QChar_exposer.def( bp::init< uchar, uchar >(( bp::arg("c"), bp::arg("r") )) );
        QChar_exposer.def( bp::init< ushort >(( bp::arg("rc") )) );
        QChar_exposer.def( bp::init< short int >(( bp::arg("rc") )) );
        QChar_exposer.def( bp::init< uint >(( bp::arg("rc") )) );
        QChar_exposer.def( bp::init< int >(( bp::arg("rc") )) );
        { //::QChar::cell
        
            typedef ::uchar ( ::QChar::*cell_function_type )(  ) const;
            cell_function_type cell_function_value( &::QChar::cell );
            
            QChar_exposer.def( 
                "cell"
                , cell_function_value );
        
        }
        { //::QChar::combiningClass
        
            typedef unsigned char ( ::QChar::*combiningClass_function_type )(  ) const;
            combiningClass_function_type combiningClass_function_value( &::QChar::combiningClass );
            
            QChar_exposer.def( 
                "combiningClass"
                , combiningClass_function_value );
        
        }
        { //::QChar::decomposition
        
            typedef ::QString ( ::QChar::*decomposition_function_type )(  ) const;
            decomposition_function_type decomposition_function_value( &::QChar::decomposition );
            
            QChar_exposer.def( 
                "decomposition"
                , decomposition_function_value );
        
        }
        { //::QChar::digitValue
        
            typedef int ( ::QChar::*digitValue_function_type )(  ) const;
            digitValue_function_type digitValue_function_value( &::QChar::digitValue );
            
            QChar_exposer.def( 
                "digitValue"
                , digitValue_function_value );
        
        }
        { //::QChar::fromAscii
        
            typedef ::QChar ( *fromAscii_function_type )( char );
            fromAscii_function_type fromAscii_function_value( &::QChar::fromAscii );
            
            QChar_exposer.def( 
                "fromAscii"
                , fromAscii_function_value
                , ( bp::arg("c") ) );
        
        }
        { //::QChar::fromLatin1
        
            typedef ::QChar ( *fromLatin1_function_type )( char );
            fromLatin1_function_type fromLatin1_function_value( &::QChar::fromLatin1 );
            
            QChar_exposer.def( 
                "fromLatin1"
                , fromLatin1_function_value
                , ( bp::arg("c") ) );
        
        }
        { //::QChar::hasMirrored
        
            typedef bool ( ::QChar::*hasMirrored_function_type )(  ) const;
            hasMirrored_function_type hasMirrored_function_value( &::QChar::hasMirrored );
            
            QChar_exposer.def( 
                "hasMirrored"
                , hasMirrored_function_value );
        
        }
        { //::QChar::isDigit
        
            typedef bool ( ::QChar::*isDigit_function_type )(  ) const;
            isDigit_function_type isDigit_function_value( &::QChar::isDigit );
            
            QChar_exposer.def( 
                "isDigit"
                , isDigit_function_value );
        
        }
        { //::QChar::isLetter
        
            typedef bool ( ::QChar::*isLetter_function_type )(  ) const;
            isLetter_function_type isLetter_function_value( &::QChar::isLetter );
            
            QChar_exposer.def( 
                "isLetter"
                , isLetter_function_value );
        
        }
        { //::QChar::isLetterOrNumber
        
            typedef bool ( ::QChar::*isLetterOrNumber_function_type )(  ) const;
            isLetterOrNumber_function_type isLetterOrNumber_function_value( &::QChar::isLetterOrNumber );
            
            QChar_exposer.def( 
                "isLetterOrNumber"
                , isLetterOrNumber_function_value );
        
        }
        { //::QChar::isLower
        
            typedef bool ( ::QChar::*isLower_function_type )(  ) const;
            isLower_function_type isLower_function_value( &::QChar::isLower );
            
            QChar_exposer.def( 
                "isLower"
                , isLower_function_value );
        
        }
        { //::QChar::isMark
        
            typedef bool ( ::QChar::*isMark_function_type )(  ) const;
            isMark_function_type isMark_function_value( &::QChar::isMark );
            
            QChar_exposer.def( 
                "isMark"
                , isMark_function_value );
        
        }
        { //::QChar::isNull
        
            typedef bool ( ::QChar::*isNull_function_type )(  ) const;
            isNull_function_type isNull_function_value( &::QChar::isNull );
            
            QChar_exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        { //::QChar::isNumber
        
            typedef bool ( ::QChar::*isNumber_function_type )(  ) const;
            isNumber_function_type isNumber_function_value( &::QChar::isNumber );
            
            QChar_exposer.def( 
                "isNumber"
                , isNumber_function_value );
        
        }
        { //::QChar::isPrint
        
            typedef bool ( ::QChar::*isPrint_function_type )(  ) const;
            isPrint_function_type isPrint_function_value( &::QChar::isPrint );
            
            QChar_exposer.def( 
                "isPrint"
                , isPrint_function_value );
        
        }
        { //::QChar::isPunct
        
            typedef bool ( ::QChar::*isPunct_function_type )(  ) const;
            isPunct_function_type isPunct_function_value( &::QChar::isPunct );
            
            QChar_exposer.def( 
                "isPunct"
                , isPunct_function_value );
        
        }
        { //::QChar::isSpace
        
            typedef bool ( ::QChar::*isSpace_function_type )(  ) const;
            isSpace_function_type isSpace_function_value( &::QChar::isSpace );
            
            QChar_exposer.def( 
                "isSpace"
                , isSpace_function_value );
        
        }
        { //::QChar::isSymbol
        
            typedef bool ( ::QChar::*isSymbol_function_type )(  ) const;
            isSymbol_function_type isSymbol_function_value( &::QChar::isSymbol );
            
            QChar_exposer.def( 
                "isSymbol"
                , isSymbol_function_value );
        
        }
        { //::QChar::isUpper
        
            typedef bool ( ::QChar::*isUpper_function_type )(  ) const;
            isUpper_function_type isUpper_function_value( &::QChar::isUpper );
            
            QChar_exposer.def( 
                "isUpper"
                , isUpper_function_value );
        
        }
        { //::QChar::mirroredChar
        
            typedef ::QChar ( ::QChar::*mirroredChar_function_type )(  ) const;
            mirroredChar_function_type mirroredChar_function_value( &::QChar::mirroredChar );
            
            QChar_exposer.def( 
                "mirroredChar"
                , mirroredChar_function_value );
        
        }
        { //::QChar::row
        
            typedef ::uchar ( ::QChar::*row_function_type )(  ) const;
            row_function_type row_function_value( &::QChar::row );
            
            QChar_exposer.def( 
                "row"
                , row_function_value );
        
        }
        { //::QChar::setCell
        
            typedef void ( ::QChar::*setCell_function_type )( ::uchar ) ;
            setCell_function_type setCell_function_value( &::QChar::setCell );
            
            QChar_exposer.def( 
                "setCell"
                , setCell_function_value
                , ( bp::arg("cell") ) );
        
        }
        { //::QChar::setRow
        
            typedef void ( ::QChar::*setRow_function_type )( ::uchar ) ;
            setRow_function_type setRow_function_value( &::QChar::setRow );
            
            QChar_exposer.def( 
                "setRow"
                , setRow_function_value
                , ( bp::arg("row") ) );
        
        }
        { //::QChar::toLower
        
            typedef ::QChar ( ::QChar::*toLower_function_type )(  ) const;
            toLower_function_type toLower_function_value( &::QChar::toLower );
            
            QChar_exposer.def( 
                "toLower"
                , toLower_function_value );
        
        }
        { //::QChar::toUpper
        
            typedef ::QChar ( ::QChar::*toUpper_function_type )(  ) const;
            toUpper_function_type toUpper_function_value( &::QChar::toUpper );
            
            QChar_exposer.def( 
                "toUpper"
                , toUpper_function_value );
        
        }
        QChar_exposer.staticmethod( "fromAscii" );
        QChar_exposer.staticmethod( "fromLatin1" );
        QChar_exposer.def("__str__", &convert );
        QChar_exposer.def("__repr__", &convert );
        QChar_exposer.def( "__copy__", &__copy__);
        QChar_exposer.def( "__deepcopy__", &__copy__);
        QChar_exposer.def( "clone", &__copy__);
        QChar_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::QChar >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        QChar_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::QChar >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
    }

}
