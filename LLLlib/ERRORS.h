#include <utils/CxErrors.hpp>
#include <utils/CxBadValueError.hpp>
#include <utils/CxIndextoSmallError.hpp>
#include <utils/CxErrors.hpp>
#include <utils/CxNullPointerError.hpp>
#include <utils/CxFileNotFoundError.hpp>
#include <utils/CxOutOfBoundsError.hpp>
#include <utils/CxIndextoLargeError.hpp>

#define ERRORTHROW(msg) throw CxErrors(msg,__FILE__,__LINE__);
#define OUTOFBOUNDSTHROW(msg) throw CxOutOfBoundsError(msg); 
#define FILENOTFOUNDTHROW(fileName) throw CxFileNotFoundError(fileName,__FILE__,__LINE__);

#define THROWERROR throw CxErrors(__FILE__,__LINE__);

