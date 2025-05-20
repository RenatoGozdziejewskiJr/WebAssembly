#ifndef MODELS_H
#define MODELS_H

#include <iostream>
#include <vector>

#include <boost/numeric/odeint.hpp>
#include <boost/thread/thread.hpp>

using namespace boost::numeric::odeint;

using state_type = std::vector<double>;

//observer
struct IntegrationObserver
{
	std::vector< state_type >& m_states;
	std::vector< double >& m_times;

	IntegrationObserver(std::vector< state_type >& states, std::vector< double >& times)
		: m_states(states), m_times(times) { }

	void operator()(const state_type& S, double t)
	{
		m_states.push_back(S);
		m_times.push_back(t);

		static unsigned long i = 0;
		std::cout << "i: " << i << " time: " << t << " : " << " y: " << S[0] << " y': " << S[1] << std::endl;
		++i;
	}
};

class IModel {
public:
	virtual ~IModel() = default;

	//The functor which the model should be implemented.
	virtual void operator() (const state_type& S, state_type& dSdt, const double t) = 0;

	virtual void setInitialState(const state_type& Si) { m_Si = Si; }
	virtual state_type getInitialState() const { return m_Si; }

private:	
	state_type m_Si; //initial State.
}; 
  
class SecondOrderModel: public IModel
{
public:
	using IModel::IModel; //Inherit all constructors from base
	
	// Inherited via Model
	void operator()(const state_type& S, state_type& dSdt, const double t) override {
		// Kp = process gain
		// taus = second order time constant
		// zeta = damping factor
		// ts ^ 2 dy2 / dt2 + 2 zeta taus dydt + y = Kp u(t - thetap)
		// taus^2*y" + 2*zeta*taus*y' + y = Kp*u(t - thetap)
		// y" = (-2*zeta*taus*y' - y + Kp*u(t - thetap)) / taus^2
		// y' = x
		// x' = (-2*zeta*taus*x - y + Kp*u(t - thetap)) / taus^2
		double u = 0.0;
		if (t > 1.0) { u = 1; }
		dSdt[0] = S[1];
		dSdt[1] = (-2.0 * m_zeta * m_taus * S[1] - S[0] + m_kp * u) / pow(m_taus, 2);
	}

	double getZeta() const { return m_zeta; }
	void setZeta(double val) { if (m_zeta != val) { m_zeta = val; } }
	 
private:
	//default process model
	//zeta = 2.0 -> overdamped step response
	//zeta = 1.0 -> critically damped step response
	//zeta = 0.5 -> underdamped step response

	double m_kp = 2.0;
	double m_taus = 1.0;
	double m_zeta = 0.5;
 
};

//Classe template que aceita somente tipos derivados de Model
template <typename T, typename = typename std::enable_if<std::is_base_of_v<IModel, T>>>
class Solver {
public:
	Solver() = default;
	T& model() { return m_model; }
	
	void setInitialTime(double t) { m_t_ini = t; }
	double getInitialTime() const { return m_t_ini; }

	void setFinalTime(double t) { m_t_end = t; }
	double getFinalTime() const { return m_t_end; }

	void setStepSize(double dt) { m_dt = dt; }
	double getStepSize() const { return m_dt; }

	std::vector<state_type> getStates() const { return m_S_vec; }
	std::vector<double> getTimes() const { return m_times_vec; }

	void clearStates() {

		m_S_vec.clear();

		m_times_vec.clear();
	}
	 
	size_t convenienceIntegration() {

		clearStates();

		state_type Si = m_model.getInitialState();

		size_t steps = integrate(m_model, Si, m_t_ini, m_t_end, m_dt, IntegrationObserver(m_S_vec, m_times_vec));

		return steps;
	}

	size_t stepperIntegrationRKDopri5() {
		
		clearStates();

		runge_kutta_dopri5< state_type > stepper;
		
		double t = m_t_ini;    // Tempo inicial

		state_type S = m_model.getInitialState();

		addStateTime(S, t);

		while (t <= m_t_end) {

			stepper.do_step(m_model, S, t, m_dt);

			t += m_dt;			

			addStateTime(S, t);

		}

		return m_S_vec.size();
	}

	size_t controlledStepperIntegration() {
		
		clearStates();

		using stepper_type = controlled_runge_kutta<runge_kutta_dopri5<state_type>>;

		stepper_type stepper;
				
		double t = m_t_ini;    // Tempo inicial
		state_type S = m_model.getInitialState();

		addStateTime(S, t);

		while (t <= m_t_end) {			
			auto res = stepper.try_step(m_model, S, t, m_dt);
			if (res == success) {
				addStateTime(S, t);
			}
		}

		return m_S_vec.size()-1; //doesnt consider the initial condition as step
	}

private:

	void addStateTime(const state_type& S, double t) {
		m_S_vec.push_back(S);

		m_times_vec.push_back(t);

		std::cout << "time: " << t << " : " << "y: " << S[0] << " y': " << S[1] << " dt: " << m_dt << std::endl;
	}

	T m_model;

	std::vector<state_type> m_S_vec; 
	std::vector<double> m_times_vec;
	double m_t_ini; // Tempo inicial
	double m_t_end;	// Tempo final
	double m_dt;	// Tamanho inicial do passo	 

};

#endif // !MODELS_H