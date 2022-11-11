%{
#include "IMP/bff/DecayRoutines.h"
%}

%apply
(double* INPLACE_ARRAY1, int DIM1) {
    (double* fit, int n_fit),
    (double *model, int n_model),
    (double* decay, int n_decay),
    (double *time_axis, int n_time_axis),
    (double *lifetime_spectrum, int n_lifetime_spectrum),
    (double *data, int n_data),
    (double* w_sq, int n_wsq),
    (double* x, int n_x),
    (double* irf, int n_irf),
    (double* irf_shift, int n_irf_shift)
}

%ignore IMP::bff::decay_rescale;
%ignore IMP::bff::decay_rescale_w;
%ignore IMP::bff::decay_rescale_w_bg;
%ignore IMP::bff::decay_fconv;
%ignore IMP::bff::decay_fconv_avx;
%ignore IMP::bff::decay_fconv_per;
%ignore IMP::bff::decay_fconv_per_avx;
%ignore IMP::bff::decay_fconv_per_cs;
%ignore IMP::bff::decay_fconv_ref;
%ignore IMP::bff::decay_sconv;
%ignore IMP::bff::decay_shift_lamp;
%inline %{
    
    double decay_rescale(
        double* fit, int n_fit,
        double* decay, int n_decay,
        int start = 0, int stop = -1
    ){
        stop = IMP::bff::mod_p(stop, n_decay + 1);
        start = IMP::bff::mod_p(start, n_decay + 1);

        double scale = 0.0;
        IMP::bff::decay_rescale(fit, decay, &scale, start, stop);
        return scale;
    }

    double decay_rescale_w(
            double* fit, int n_fit,
            double* decay, int n_decay,
            double* w_sq, int n_wsq,
            int start = 0, int stop = -1
    ){
        stop = IMP::bff::mod_p(stop, n_decay + 1);
        start = IMP::bff::mod_p(start, n_decay + 1);

        double scale = 0.0;
        IMP::bff::decay_rescale_w(fit, decay, w_sq, &scale, start, stop);
        return scale;
    }

    double decay_rescale_w_bg(
            double* fit, int n_fit,
            double* decay, int n_decay,
            double* w_sq, int n_wsq,
            double bg = 0.0,
            int start = 0, int stop = -1
    ){
        stop = IMP::bff::mod_p(stop, n_decay + 1);
        start = IMP::bff::mod_p(start, n_decay + 1);

        double scale = 0.0;
        IMP::bff::decay_rescale_w_bg(fit, decay, w_sq, bg, &scale, start, stop);
        return scale;
    }

    void decay_fconv(
            double* fit, int n_fit,
            double* irf, int n_irf,
            double* x, int n_x,
            int start = 0, int stop = -1,
            double dt = 1.0
    ){
        stop = IMP::bff::mod_p(stop, n_fit + 1);
        start = IMP::bff::mod_p(start, n_fit + 1);
        IMP::bff::decay_fconv(fit, x, irf, n_x / 2, start, stop, dt);
    }

    void decay_fconv_avx(
            double* fit, int n_fit,
            double* irf, int n_irf,
            double* x, int n_x,
            int start = 0, int stop = -1,
            double dt = 1.0
    ){
        stop = IMP::bff::mod_p(stop, n_fit + 1);
        start = IMP::bff::mod_p(start, n_fit + 1);
        IMP::bff::decay_fconv_avx(fit, x, irf, n_x / 2, start, stop, dt);
    }

    void decay_fconv_per(
            double* fit, int n_fit,
            double* irf, int n_irf,
            double* x, int n_x,
            double period,
            int start = 0, int stop = -1,
            double dt = 1.0
    ){
        stop = IMP::bff::mod_p(stop, n_fit + 1);
        start = IMP::bff::mod_p(start, n_fit + 1);
        IMP::bff::decay_fconv_per(fit, x, irf, n_x / 2, start, stop, n_fit, period, dt);
    }

    void decay_fconv_per_avx(
            double* fit, int n_fit,
            double* irf, int n_irf,
            double* x, int n_x,
            double period,
            int start = 0, int stop = -1,
            double dt = 1.0
    ){
        stop = IMP::bff::mod_p(stop, n_fit + 1);
        start = IMP::bff::mod_p(start, n_fit + 1);
        IMP::bff::decay_fconv_per_avx(fit, x, irf, n_x / 2, start, stop, n_fit, period, dt);
    }

    void decay_fconv_per_cs(
            double* fit, int n_fit,
            double* irf, int n_irf,
            double* x, int n_x,
            double period,
            int conv_stop = -1,
            int stop = -1,
            double dt = 1.0
    ){
        stop = IMP::bff::mod_p(stop, n_fit + 1);
        IMP::bff::decay_fconv_per_cs(fit, x, irf, n_x / 2, stop, n_fit, period, conv_stop, dt);
    }

    void decay_fconv_ref(
            double* fit, int n_fit,
            double* irf, int n_irf,
            double* x, int n_x,
            double tauref,
            int start = 0, int stop = -1,
            double dt = 1.0
    ){
        stop = IMP::bff::mod_p(stop, n_fit + 1);
        start = IMP::bff::mod_p(start, n_fit + 1);

        IMP::bff::decay_fconv_ref(fit, x, irf, n_x / 2, start, stop, tauref);
    }

    void decay_sconv(
            double* fit, int n_fit,
            double* irf, int n_irf,
            double* model, int n_model,
            int start = 0, int stop = -1
    ){
        stop = IMP::bff::mod_p(stop, n_fit + 1);
        start = IMP::bff::mod_p(start, n_fit + 1);
        IMP::bff::decay_sconv(fit, model, irf, start, stop);
    }

    void decay_shift_lamp(
            double* irf, int n_irf,
            double* irf_shift, int n_irf_shift,
            double ts
    ){
        IMP::bff::decay_shift_lamp(irf_shift, irf, ts, n_irf);
    }

%}

//%include "IMP/bff/DecayRoutines.h"



