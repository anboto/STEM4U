#ifndef _STEM4U_Iir_fir_h_
#define _STEM4U_Iir_fir_h_

namespace Upp {
	
template <typename T>
class Linear_IIR_Filter {
public:
	Linear_IIR_Filter() {}
    Linear_IIR_Filter(T alpha) : alpha(alpha) {}

    T Filter(T input) {		// y[n] = αx[n] + (1−α)y[n−1]
        if (IsNull(prevOutput))
            prevOutput = input;
        return prevOutput = alpha*input + (1 - alpha)*prevOutput;
    }
	void LowPassCoefficients(T cutoffFrequency /* Hz */, T samplingPeriod /* sec */) {
		T omega_c = 2*M_PI*cutoffFrequency*samplingPeriod;
		alpha = omega_c/(1 + omega_c);
	}
	void HighPassCoefficients(T cutoffFrequency /* Hz */, T samplingPeriod /* sec */) {
		T omega_c = M_PI*cutoffFrequency*samplingPeriod;
		T tan_omega_c = tan(omega_c);
    	alpha = (1 - tan_omega_c)/(1 + tan_omega_c);
	}
private:
    T alpha;	// α = Δt/(RC + Δt)
    T prevOutput = Null;
};

template <typename T>
class BiquadFilter {
public:
	BiquadFilter() {}
    BiquadFilter(T b0, T b1, T b2, T a1, T a2) : b0(b0), b1(b1), b2(b2), a1(a1), a2(a2) {}
	
    T Filter(T input) {
        if (IsNull(x1))
            x1 = x2 = y1 = y2 = input;
        
        T output = b0*input + b1*x1 + b2*x2 - a1*y1 - a2*y2;

        // Update previous samples
        x2 = x1;
        x1 = input;
        y2 = y1;
        y1 = output;

        return output;
    }
    // q ia a quality factor. Good values may be between 1/2 to 1/sqrt(2) and 1
    // As q increases, the roll-off becomes steeper, and there may be a peak at the cutoff frequency
    void LowPassCoefficients(T cutoffFrequency, T samplingPeriod, T q) {
        T omega_c = 2*M_PI*cutoffFrequency*samplingPeriod;
        T alpha = sin(omega_c)/2/q;
        T cos_omega_c = cos(omega_c);

        T a0 = 1 + alpha;
        b0 = (1 - cos_omega_c)/2/a0;
        b1 = (1 - cos_omega_c)/a0;
        b2 = (1 - cos_omega_c)/2/a0;
        a1 = -2*cos_omega_c/a0;
        a2 = (1 - alpha)/a0;
    }
    void HighPassCoefficients(T cutoffFrequency, T samplingPeriod, T q) {
        T omega_c = 2*M_PI*cutoffFrequency*samplingPeriod;
        T alpha = sin(omega_c)/2/q;
        T cos_omega_c = cos(omega_c);

        T a0 = 1 + alpha;
        b0 = (1 + cos_omega_c)/2/a0;
        b1 = -(1 + cos_omega_c)/a0;
        b2 = (1 + cos_omega_c)/2/a0;
        a1 = -2*cos_omega_c/a0;
        a2 = (1 - alpha)/a0;
    }
    // For band-pass filters, q is related to the bandwidth of the passband
    // Higher q narrows the passband
    // q = 1/sqrt(2) could be a good intermediate value
    void BandPassCoefficients(T cutoffFrequency, T samplingPeriod, T q) {
        T omega_c = 2*M_PI*cutoffFrequency*samplingPeriod;
        T alpha = sin(omega_c)/2/q;
        T cos_omega_c = cos(omega_c);

        // Coefficients for a low-pass filter
        T a0 = 1 + alpha;
        b0 = alpha/a0;
        b1 = 0;
        b2 = -b0;
        a1 = -2*cos_omega_c/a0;
        a2 = (1 - alpha)/a0;
    }
    
private:
    T b0, b1, b2;
    T a1, a2;
    T x1 = Null, x2 = Null; // Previous inputs
    T y1 = Null, y2 = Null; // Previous outputs
};

template <typename T>
class FIRFilter {
public:
	FIRFilter() {}
	FIRFilter(int windowWidth) {
		kernel.resize(windowWidth, 1./windowWidth);
	}
	void LowPassCoefficients(T cutoffFrequency, T samplingPeriod) {
		int windowWidth = int(1/cutoffFrequency/samplingPeriod) + 1;	
		kernel = Eigen::Matrix<T, Dynamic, 1>::Constant(windowWidth, 1./windowWidth);
	}
	T Filter(T input) {
		if (buffer.size() == 0)
			buffer = Eigen::Matrix<T, Dynamic, 1>::Constant(kernel.size(), input);
		Eigen::Index len = buffer.size()-1;
		buffer.tail(len) = buffer.head(len).eval();
		buffer(0) = input;
		return buffer.cwiseProduct(kernel).sum();
	}
	FIRFilter& operator=(const FIRFilter& f) {
		buffer = f.buffer;
		kernel = f.kernel;	
		return *this;
	}
	
private:
	Eigen::Matrix<T, Dynamic, 1> buffer, kernel;
};

}

#endif
