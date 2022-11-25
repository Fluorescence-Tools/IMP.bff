#include <IMP/bff/DecayLinearization.h>

IMPBFF_BEGIN_NAMESPACE

void DecayLinearization::set_linearization_table(DecayCurve* v) {
    set_data(v);
}

DecayCurve* DecayLinearization::get_linearization_table() {
    return get_data();
}

void DecayLinearization::add(DecayCurve *out) {
    if ((out != nullptr) && is_active()) {
        DecayCurve *lt = get_linearization_table();
        lt->resize(out->size(), 1.0);
        size_t start, stop;
        start = std::min(get_start(out), out->size());
        stop = std::min(get_stop(out), out->size());
        for (size_t i = start; i < stop; i++) {
            out->y[i] *= lt->y[i];
        }
    }
}

DecayLinearization::DecayLinearization(
        DecayCurve* linearization_table,
        int start, int stop,
        bool active,
        int n_window
) : DecayModifier(nullptr, start, stop, active){
#if IMPBFF_VERBOSE
    std::clog << "DecayLinearization::DecayLinearization(double*, int, bool)" << std::endl;
#endif
    linearization_table_.set_x(linearization_table->get_x());
    linearization_table_.set_y(linearization_table->get_y());
    linearization_table_.apply_simple_moving_average(
            get_start(linearization_table),
            get_stop(linearization_table),
            n_window,
            true
    );
    set_linearization_table(&linearization_table_);
    set_start(start);
    set_stop(stop);
}

IMPBFF_END_NAMESPACE
