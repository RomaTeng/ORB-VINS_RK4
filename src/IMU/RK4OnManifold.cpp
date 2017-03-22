#include "RK4OnManifold.h"

namespace ORB_SLAM2
{

RK4Preint_State::RK4Preint_State( )
{
	Reset( );
}

RK4Preint_State::RK4Preint_State( const RK4Preint_State& other )
	: _J1( other._J1 )
	, _J2( other._J2 )
	, _J3( other._J3 )
	, _j1_w( other._j1_w )
	, _j2_w( other._j2_w )
	, _j2_a( other._j2_a )
	, _j3_w( other._j3_w )
	, _j3_a( other._j3_a )
{
}

RK4Preint_State&
RK4Preint_State::operator+=( const K_State & other )
{
	_J1 *= Sophus::SO3::expMatrix( other._kj_dJ1 );
    _J1 =  Sophus::SO3(_J1).matrix();
	_J2 += other._kj_dJ2;
	_J3 += other._kj_dJ3;

	_j1_w += other._kj_dJ1_w;
	_j2_a += other._kj_dJ2_a;
	_j2_w += other._kj_dJ2_w;
	_j3_a += other._kj_dJ3_a;
	_j3_w += other._kj_dJ3_w;

	return *this;
}

void
RK4Preint_State::Write(std::ostream& os) const
{
	os  << "J1: " << _J1 << std::endl
		<< "J2: " << _J2 << std::endl
		<< "J3: " << _J3 << std::endl
		<< "j1_w: " << _j1_w  << std::endl
		<< "j2_w: " << _j2_w  << std::endl
		<< "j2_a: " << _j2_a  << std::endl
		<< "j3_w: " << _j3_w  << std::endl
		<< "j3_a: " << _j3_a  << std::endl;
}

RK4Preint_Cov::RK4Preint_Cov()
{
	Reset( );
}

RK4Preint_Cov::RK4Preint_Cov( const RK4Preint_Cov& other )
	: _P( other._P )
{
}

RK4Preint_Cov&
RK4Preint_Cov::operator+=( const K_Cov & other )
{
	_P += other._kj_P;
	return *this;
}

void
RK4Preint_Cov::Write(std::ostream& os) const
{
	os  << "P: " << _P << std::endl;
}

RK4Preint::RK4Preint( )
	: _J_State( )
	, _J_Cov( )
	, _prev_dt( 0 )
	, _isFirst( true )
{
}

RK4Preint::RK4Preint( const RK4Preint& other )
	: _J_State( other._J_State )
	, _J_Cov( other._J_Cov )
	, _prev_dt( other._prev_dt )
	, _isFirst( other._isFirst )
{
}

void
RK4Preint::Reset( )
{
	_isFirst = true;
	_J_State.Reset( );
	_J_Cov.Reset( );
}

void
RK4Preint::IntegrateOnce( const Vector3d& omega, const Vector3d& acc, const double dt, const Matrix9d& noiseCov )
{
	double h = dt;

	if ( _isFirst )
	{
		_isFirst = false;

		{
			K_State k( omega, acc, _J_State, h );
			_J_State += k;
		}
		{
			K_Cov K( omega, acc, noiseCov, _J_Cov, h );
			_J_Cov += K;
		}

		UpdateMemory( omega, acc, dt );
	}
	else
	{
		Vector3d	omega_predicted1, omega_predicted2, acc_predicted1, acc_predicted2;
		if ( _prev_dt < 0.000001 )
		{
			omega_predicted1 =  omega;
			omega_predicted2 =  omega;
			acc_predicted1   =  acc;
			acc_predicted2   =  acc;
		}
		else
		{
			omega_predicted1 =  ( omega-_prev_omega ) * dt / _prev_dt * 0.5 + omega;
		    omega_predicted2 =  ( omega-_prev_omega ) * dt / _prev_dt *   1 + omega;
		    acc_predicted1   =  ( acc  -_prev_acc)    * dt / _prev_dt * 0.5 + acc;
		    acc_predicted2   =  ( acc  -_prev_acc)    * dt / _prev_dt *   1 + acc;
		}
		// RK4 for State and Jacobians
		{

			K_State k1( omega, acc, _J_State, h );
			//RK4Preint y1 = _J_State + k1 * 0.5;

			K_State k2( omega_predicted1, acc_predicted1, _J_State + k1 * 0.5, h );
			//RK4Preint y2 = _J_State + k2 * 0.5;

			K_State k3( omega_predicted1, acc_predicted1, _J_State + k2 * 0.5, h );
			//RK4Preint y3 = _J_State + k3;

			K_State k4( omega_predicted2, acc_predicted2, _J_State + k3, h );

			// ( ( k1 + 2 * k2 + 2 * k3 + k4 ) / 6 )
			k1 *= 1./6;
			k2 *= 2./6;
			k3 *= 2./6;
			k4 *= 1./6;

			k1 += ( k2 += ( k3 += k4 ) );

			_J_State += ( k1 );
		}

		// RK4 for Covariance
		{
			// covariance
			K_Cov K1( omega, acc, noiseCov, _J_Cov, h );
			//Matrix9d P1 = _J_Cov + K1._kj_P * 0.5;

			K_Cov K2( omega_predicted1, acc_predicted1, noiseCov, _J_Cov + K1 * 0.5, h );
			//Matrix9d P2 = _P + K2._kj_P * 0.5;

			K_Cov K3( omega_predicted1, acc_predicted1, noiseCov, _J_Cov + K2 * 0.5, h );
			//Matrix9d P3 = _P + K3._kj_P ;

			K_Cov K4( omega_predicted2, acc_predicted2, noiseCov, _J_Cov + K3, h );

			K1 *= 1./6;
			K2 *= 2./6;
			K3 *= 2./6;
			K4 *= 1./6;

			K1 += ( K2 += ( K3 += K4 ) );

			_J_Cov += ( K1 );
		}

		UpdateMemory( omega, acc, dt );

	}
}

void
RK4Preint::Write(std::ostream& os) const
{
	_J_State.Write( os );
	_J_Cov.Write( os );
}

K_State::K_State( const Vector3d& omega, const Vector3d& acc, const RK4Preint_State& JState, double h )
{
	_kj_dJ1 = omega;
	_kj_dJ2 = JState._J1 * acc;
	_kj_dJ3 = JState._J2;

	_kj_dJ1_w = - Sophus::SO3::hat( omega ) * JState._j1_w - Matrix3d::Identity( );
	_kj_dJ2_w = - JState._J1 * Sophus::SO3::hat( acc ) * JState._j1_w;
	_kj_dJ2_a = - JState._J1;
	_kj_dJ3_w = JState._j2_w;
	_kj_dJ3_a = JState._j2_a;

	(*this) *= h;
}

K_State&
K_State::operator*=( const double h )
{
	_kj_dJ1 *= h;
	_kj_dJ2 *= h;
	_kj_dJ3 *= h;

	_kj_dJ1_w *= h;
	_kj_dJ2_w *= h;
	_kj_dJ2_a *= h;
	_kj_dJ3_w *= h;
	_kj_dJ3_a *= h;

	return *this;
}

K_State&
K_State::operator+=( const K_State & other )
{
	_kj_dJ1 += other._kj_dJ1;
	_kj_dJ2 += other._kj_dJ2;
	_kj_dJ3 += other._kj_dJ3;

	_kj_dJ1_w += other._kj_dJ1_w;
	_kj_dJ2_a += other._kj_dJ2_a;
	_kj_dJ2_w += other._kj_dJ2_w;
	_kj_dJ3_a += other._kj_dJ3_a;
	_kj_dJ3_w += other._kj_dJ3_w;

	return *this;
}

void
K_State::Write(std::ostream& os) const
{
	os  << "kj_dJ1: " << _kj_dJ1 << std::endl
		<< "kj_dJ2: " << _kj_dJ2 << std::endl
		<< "kj_dJ1_w: " << _kj_dJ1_w  << std::endl
		<< "kj_dJ2_w: " << _kj_dJ2_w  << std::endl
		<< "kj_dJ2_a: " << _kj_dJ2_a  << std::endl
		<< "kj_dJ3_w: " << _kj_dJ3_w  << std::endl
		<< "kj_dJ3_a: " << _kj_dJ3_a  << std::endl;
}

K_Cov::K_Cov( const Vector3d& omega, const Vector3d& acc, const Matrix9d& noiseCov, const RK4Preint_Cov& JCov, double h )
{
	Matrix9d A;
	A.setZero();

	Matrix3d tempWHat = - Sophus::SO3::hat( omega );

	A.block<3,3>(0,0) = tempWHat;
	A.block<3,3>(3,0) = - Sophus::SO3::hat( acc );
	A.block<3,3>(3,3) = tempWHat;
	A.block<3,3>(6,3) = Matrix3d::Identity( );
	A.block<3,3>(6,6) = tempWHat;

	const Matrix9d& P0 = JCov._P;
	_kj_P = ( A * P0 + P0 * A.transpose() + noiseCov ) * h;
}

K_Cov&
K_Cov::operator*=( const double h )
{
	_kj_P *= h;
	return *this;
}

K_Cov&
K_Cov::operator+=( const K_Cov& other )
{
	_kj_P += other._kj_P;
	return *this;
}

void
K_Cov::Write(std::ostream& os) const
{
	os << "kj_P: " << _kj_P << std::endl;
}

RK4Preint_State
operator+( const RK4Preint_State& Y, const K_State& K )
{
	RK4Preint_State yOut = Y;
	yOut += K;
	return yOut;
}

RK4Preint_Cov
operator+( const RK4Preint_Cov& P, const K_Cov& KP )
{
	RK4Preint_Cov pOut = P;
	pOut += KP;
	return pOut;
}

K_State
operator*( const K_State& K, const double h )
{
	K_State kOut = K;
	kOut *= h;
	return kOut;
}

K_Cov
operator*( const K_Cov& P, const double h )
{
	K_Cov pOut = P;
	pOut *= h;
	return pOut;
}

}
