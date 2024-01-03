#ifndef PRAGMA_H
#define PRAGMA_H

#define LSCALAR float
#define epsDir 1e-3f
#define epsDist 1e-3f

#if defined(_MSC_VER) and !defined(__CUDACC__)
#define _USE_MATH_DEFINES // for C++
#define DLL_EXPORT __declspec(dllexport)
#include <cmath>
#else
#define DLL_EXPORT
#endif

//cuda macro
#ifdef __CUDACC__
extern void initializeGPU();
extern void finalizeGPU();
#define DEVICE_HOST __device__ __host__
#define FINALIZE_GPU finalizeGPU();
#else
#define DEVICE_HOST
#define FINALIZE_GPU
#endif

//openmp macro prepare
#ifdef _MSC_VER
#define STRINGIFY_OMP(X) X
#define PRAGMA __pragma
#else
#define STRINGIFY_OMP(X) #X
#define PRAGMA _Pragma
#endif
#ifndef NO_OPENMP
//openmp convenient functions
#define OMP_PARALLEL_FOR_ PRAGMA(STRINGIFY_OMP(omp parallel for schedule(static)))
#define OMP_PARALLEL_FOR_I(...) PRAGMA(STRINGIFY_OMP(omp parallel for schedule(static) __VA_ARGS__))
#define OMP_PARALLEL_FOR_X(X) PRAGMA(STRINGIFY_OMP(omp parallel for num_threads(X) schedule(static)))
#define OMP_PARALLEL_FOR_XI(X,...) PRAGMA(STRINGIFY_OMP(omp parallel for num_threads(X) schedule(static) __VA_ARGS__))
#define OMP_ADD(...) reduction(+: __VA_ARGS__)
#define OMP_PRI(...) private(__VA_ARGS__)
#define OMP_FPRI(...) firstprivate(__VA_ARGS__)
#define OMP_ATOMIC_ PRAGMA(STRINGIFY_OMP(omp atomic))
#ifdef _MSC_VER
#define OMP_ATOMIC_CAPTURE_ PRAGMA(STRINGIFY_OMP(omp critical))	// VS doesn't support capture, use critical instead
#else
#define OMP_ATOMIC_CAPTURE_ PRAGMA(STRINGIFY_OMP(omp atomic capture))
#endif
#define OMP_CRITICAL_ PRAGMA(STRINGIFY_OMP(omp critical))
#define OMP_FLUSH_(X) PRAGMA(STRINGIFY_OMP(omp flush(X)))
#else
//openmp convenient functions
#define OMP_PARALLEL_FOR_
#define OMP_PARALLEL_FOR_I(...)
#define OMP_PARALLEL_FOR_X(X)
#define OMP_PARALLEL_FOR_XI(X,...)
#define OMP_ADD(...)
#define OMP_PRI(...)
#define OMP_FPRI(...)
#define OMP_ATOMIC_
#define OMP_ATOMIC_CAPTURE_
#define OMP_CRITICAL_
#define OMP_FLUSH_(X)
#endif

//assert
#define ASSERT(var) {do{if(!(var)){FINALIZE_GPU exit(EXIT_FAILURE);}}while(0);}
#define ASSERT_MSG(var,msg) {do{if(!(var)){printf(msg);FINALIZE_GPU exit(EXIT_FAILURE);}}while(0);}
#define ASSERT_MSGV(var,msg,...) {do{if(!(var)){printf(msg,__VA_ARGS__);FINALIZE_GPU exit(EXIT_FAILURE);}}while(0);}
#define FUNCTION_NOT_IMPLEMENTED ASSERT_MSGV(false,"Function \"%s\" not implemented!",__FUNCTION__)

//basic rigid body transformation
#define TRANSI(JSS,I) JSS.template block<3,4>(0,(I)*4)
#define TRANSI6(JSS,I) JSS.template block<6,6>(0,(I)*6)
#define GETT(TID,JSS,I) Eigen::Block<Mat3Xd,3,4> TID=TRANSI(JSS,I);
#define GETTC(TID,JSS,I) Eigen::Block<const Mat3Xd,3,4> TID=TRANSI(JSS,I);
#define GETTM(TID,JSS,I) Eigen::Block<Eigen::Map<Mat3Xd >,3,4> TID=TRANSI(JSS,I);
#define GETTCM(TID,JSS,I) Eigen::Block<const Eigen::Map<Mat3Xd >,3,4> TID=TRANSI(JSS,I);
#define GETT_T(TID,JSS,I) Eigen::Block<Mat3XT,3,4> TID=TRANSI(JSS,I);
#define GETTC_T(TID,JSS,I) Eigen::Block<const Mat3XT,3,4> TID=TRANSI(JSS,I);
#define GETTM_T(TID,JSS,I) Eigen::Block<Eigen::Map<Mat3XT >,3,4> TID=TRANSI(JSS,I);
#define GETTCM_T(TID,JSS,I) Eigen::Block<const Eigen::Map<Mat3XT >,3,4> TID=TRANSI(JSS,I);
#define ROT(A) (A).template block<3,3>(0,0)
#define CTR(A) (A).template block<3,1>(0,3)
#define SCALE(A) sqrt((ROT(A).transpose()*ROT(A)).diagonal().mean())
#define ROT_NO_SCALE(A) ROT(A)/SCALE(A)
#define ROTI(A,I) (A).template block<3,3>(0,(I)*4+0)
#define CTRI(A,I) (A).template block<3,1>(0,(I)*4+3)
#define DWDLI(A,I) (A).col(((I)<<1)+0)
#define DTDLI(A,I) (A).col(((I)<<1)+1)

//6x6 spatial block
#define SPATIAL_BLK00(A) A.template block<3,3>(0,0)
#define SPATIAL_BLK10(A) A.template block<3,3>(3,0)
#define SPATIAL_BLK01(A) A.template block<3,3>(0,3)
#define SPATIAL_BLK11(A) A.template block<3,3>(3,3)

//multiply transformation
#define APPLY_TRANS(C,A,B)\
C=(ROT(A)*(B)).eval();\
CTR(C)+=CTR(A);

//invert the transformation
#define INV(IA,A)\
ROT(IA)=ROT(A).transpose().eval();\
CTR(IA)=-(ROT(IA)*CTR(A)).eval();

//degree to radius
#define D2R(X) (X*M_PI/180.0f)

#define DECL_MAT_VEC_MAP_TYPES_T \
typedef Eigen::Matrix<T,2,1> Vec2T;\
typedef Eigen::Matrix<T,3,1> Vec3T;\
typedef Eigen::Matrix<T,4,1> Vec4T;\
typedef Eigen::Matrix<int,3,1> Vec3i;\
typedef Eigen::Matrix<T,6,1> Vec6T;\
typedef Eigen::Matrix<int,9,1> Vec9i;\
typedef Eigen::Matrix<T,24,1> Vec24T;\
typedef Eigen::Matrix<T,27,1> Vec27T;\
typedef Eigen::Matrix<int,27,1> Vec27i;\
typedef Eigen::Matrix<T,-1,1> Vec;\
typedef Eigen::Matrix<int,-1,1> Veci;\
typedef Eigen::Matrix<unsigned int,-1,1> Vecui;\
\
typedef Eigen::Matrix<T,3,3> Mat3T;\
typedef Eigen::Matrix<T,3,4> Mat3X4T;\
typedef Eigen::Matrix<T,-1,6> MatX6T;\
typedef Eigen::Matrix<int,-1,8> MatX8i;\
typedef Eigen::Matrix<T,8,8> Mat8T;\
typedef Eigen::Matrix<T,24,24> Mat24T;\
typedef Eigen::Matrix<int,-1,27> MatX27i;\
typedef Eigen::Matrix<T,-1,27*9> MatXST;\
typedef Eigen::Matrix<T,-1,-1> MatT;\
\
typedef Eigen::Map<Vec2T> Vec2TM;\
typedef Eigen::Map<const Vec2T> Vec2TCM;\
typedef Eigen::Map<Vec3T> Vec3TM;\
typedef Eigen::Map<const Vec3T> Vec3TCM;\
typedef Eigen::Map<Vec24T> Vec24TM;\
typedef Eigen::Map<const Vec24T> Vec24TCM;\
typedef Eigen::Map<Vec> VecM;\
typedef Eigen::Map<const Vec> VecCM;\
typedef Eigen::Map<Vecui> VecuiM;\
typedef Eigen::Map<const Vecui> VecuiCM;\
\
typedef Eigen::Map<Mat3T> Mat3TM;\
typedef Eigen::Map<const Mat3T> Mat3TCM;\
typedef Eigen::Map<MatX8i> MatX8iM;\
typedef Eigen::Map<const MatX8i> MatX8iCM;\
typedef Eigen::Map<Mat8T> Mat8TM;\
typedef Eigen::Map<const Mat8T> Mat8TCM;\
typedef Eigen::Map<Mat24T> Mat24TM;\
typedef Eigen::Map<const Mat24T> Mat24TCM;\
typedef Eigen::Map<MatX27i> MatX27iM;\
typedef Eigen::Map<const MatX27i> MatX27iCM;\
typedef Eigen::Map<MatXST> MatXSTM;\
typedef Eigen::Map<const MatXST> MatXSTCM;\
typedef Eigen::Map<MatT> MatTM;\
typedef Eigen::Map<const MatT> MatTCM;\
\
typedef Eigen::Map<Vec3T,Eigen::Unaligned,Eigen::InnerStride<>> Vec3TSM;\
typedef Eigen::Map<const Vec3T,Eigen::Unaligned,Eigen::InnerStride<>> Vec3TSCM;\
typedef Eigen::Map<Vec27T,Eigen::Unaligned,Eigen::InnerStride<>> Vec27TSM;\
typedef Eigen::Map<const Vec27T,Eigen::Unaligned,Eigen::InnerStride<>> Vec27TSCM;\
typedef Eigen::Map<Mat3T,Eigen::Unaligned,Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>> Mat3TSM;\
typedef Eigen::Map<const Mat3T,Eigen::Unaligned,Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>> Mat3TSCM;\
\
typedef Eigen::SparseMatrix<T,0,int> SMatT;\
typedef Eigen::Triplet<T,int> STrip;\
typedef ParallelVector<STrip> STrips;

#define DECL_MAP_FUNCS  \
template <typename T2>   \
static inline Eigen::Map<T2> mapM(T2& m) {   \
  return Eigen::Map<T2>(m.data(),m.rows(),m.cols());  \
}   \
template <typename T2>   \
static inline Eigen::Map<T2> mapV(T2& m) {   \
  return Eigen::Map<T2>(m.data(),m.rows());  \
}   \
template <typename T2>   \
static inline Eigen::Map<T2> mapM(T2* m) {   \
  static const int rows=T2::RowsAtCompileTime>0?T2::RowsAtCompileTime:0;  \
  static const int cols=T2::ColsAtCompileTime>0?T2::ColsAtCompileTime:0;  \
  return m?Eigen::Map<T2>(m->data(),m->rows(),m->cols()):Eigen::Map<T2>(NULL,rows,cols);  \
}   \
template <typename T2>   \
static inline Eigen::Map<T2> mapV(T2* m) {   \
  return m?Eigen::Map<T2>(m->data(),m->rows()):Eigen::Map<T2>(NULL,0);  \
}   \
template <typename T2>   \
static inline Eigen::Map<const T2> mapCM(const T2& m) {   \
  return Eigen::Map<const T2>(m.data(),m.rows(),m.cols());  \
}   \
template <typename T2>   \
static inline Eigen::Map<const T2> mapCV(const T2& m) {   \
  return Eigen::Map<const T2>(m.data(),m.rows());  \
}   \
template <typename T2>   \
static inline Eigen::Map<const T2> mapCM(const T2* m) {   \
  static const int rows=T2::RowsAtCompileTime>0?T2::RowsAtCompileTime:0;  \
  static const int cols=T2::ColsAtCompileTime>0?T2::ColsAtCompileTime:0;  \
  return m?Eigen::Map<const T2>(m->data(),m->rows(),m->cols()):Eigen::Map<const T2>(NULL,rows,cols);  \
}   \
template <typename T2>   \
static inline Eigen::Map<const T2> mapCV(const T2* m) {   \
  return m?Eigen::Map<const T2>(m->data(),m->rows()):Eigen::Map<const T2>(NULL,0);  \
}   \
template <typename T2>   \
static inline Eigen::Map<const T2> mapM2CM(Eigen::Map<T2> m) {   \
  return Eigen::Map<const T2>(m.data(),m.rows(),m.cols());  \
}   \
template <typename T2>   \
static inline Eigen::Map<T2> mapCM2M(Eigen::Map<const T2> m) {   \
  return Eigen::Map<T2>(m.data(),m.rows(),m.cols());  \
}   \
template <typename T2>   \
static inline Eigen::Map<const T2> mapV2CV(Eigen::Map<T2> m) {   \
  return Eigen::Map<const T2>(m.data(),m.rows());  \
}   \
template <typename T2>   \
static inline Eigen::Map<T2> mapCV2V(Eigen::Map<const T2> m) {   \
  return Eigen::Map<T2>(m.data(),m.rows());  \
}

#define REUSE_MAP_FUNCS_T(T)  \
using T::mapM;\
using T::mapV;\
using T::mapCM;\
using T::mapCV;\
using T::mapM2CM;\
using T::mapCM2M;\
using T::mapV2CV;\
using T::mapCV2V;

#endif
