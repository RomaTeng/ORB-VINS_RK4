#ifndef TVISLAM_RK4ONMANIFOLD_H
#define TVISLAM_RK4ONMANIFOLD_H

#include <Eigen/Dense>

#include "imudata.h"
#include "so3.h"

namespace ORB_SLAM2
{

typedef Eigen::Matrix<double, 9, 9> Matrix9d;

struct K_State;
struct K_Cov;

struct RK4Preint_State
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	RK4Preint_State( );
	RK4Preint_State( const RK4Preint_State& other );

	inline void
	Reset( )
	{
		_J1.setIdentity();
		_J2.setZero();
		_J3.setZero();
		_j1_w.setZero();
		_j2_w.setZero();
		_j2_a.setZero();
		_j3_w.setZero();
		_j3_a.setZero();
	}

	RK4Preint_State& operator+=( const K_State& other );
	void Write(std::ostream& os) const;

	Eigen::Matrix3d _J1;
	Eigen::Vector3d _J2, _J3;
	Eigen::Matrix3d _j1_w, _j2_w, _j2_a, _j3_w, _j3_a;
};

struct RK4Preint_Cov
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	RK4Preint_Cov( );
	RK4Preint_Cov( const RK4Preint_Cov& other );

	inline void
	Reset( )
	{
		_P.setZero();
	}

	RK4Preint_Cov& operator+=( const K_Cov& other );
	void Write(std::ostream& os) const;

	Matrix9d _P;
};

struct RK4Preint
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	RK4Preint( );
	RK4Preint( const RK4Preint& other );

	void IntegrateOnce( const Vector3d& omega, const Vector3d& acc, const double dt, const Matrix9d& noise );

	void Reset( );

	void inline
	UpdateMemory( const Vector3d& omega, const Vector3d& acc, const double dt )
	{
		_prev_omega = omega;
		_prev_acc   = acc;
		_prev_dt    = dt;
	}

	void Write(std::ostream& os) const;

	RK4Preint_State				_J_State;
	RK4Preint_Cov				_J_Cov;

	Vector3d					_prev_omega, _prev_acc;
	double						_prev_dt;

	bool _isFirst;

};

struct K_State
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	K_State( const Vector3d& omega, const Vector3d& acc, const RK4Preint_State& y, double h );

    K_State& operator*=( const double h );

    K_State& operator+=( const K_State& other );

	void Write(std::ostream& os) const;

    Eigen::Vector3d				_kj_dJ1, _kj_dJ2, _kj_dJ3;
	Eigen::Matrix3d 			_kj_dJ1_w,
								_kj_dJ2_a,	_kj_dJ2_w,
								_kj_dJ3_a,	_kj_dJ3_w;

    double                      _h;
};

struct K_Cov
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	K_Cov( const Vector3d& omega, const Vector3d& acc, const Matrix9d& noise, const RK4Preint_Cov& kj_P, double h );

    K_Cov& operator*=( const double h );

    K_Cov& operator+=( const K_Cov& other );

	void Write(std::ostream& os) const;

	Matrix9d _kj_P;

};

RK4Preint_State
operator+( const RK4Preint_State& Y, const K_State& K );

RK4Preint_Cov
operator+( const RK4Preint_Cov& Y, const K_Cov& K );

K_State
operator*( const K_State& K, const double h );

K_Cov
operator*( const K_Cov& P, const double h );

inline std::ostream&
operator<<(std::ostream& s, const RK4Preint& rvpb)
{ rvpb.Write(s); return s;  }

inline std::ostream&
operator<<(std::ostream& s, const RK4Preint_State& rvpb)
{ rvpb.Write(s); return s;  }

inline std::ostream&
operator<<(std::ostream& s, const RK4Preint_Cov& rvpb)
{ rvpb.Write(s); return s;  }

inline std::ostream&
operator<<(std::ostream& s, const K_State& rvpb_K)
{ rvpb_K.Write(s); return s;  }

inline std::ostream&
operator<<(std::ostream& s, const K_Cov& rvpb_P)
{ rvpb_P.Write(s); return s;  }

} // namespace

#endif // TVISLAM_IMUPREINTEGRATOR_H
