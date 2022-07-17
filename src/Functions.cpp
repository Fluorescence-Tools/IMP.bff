#include <IMP/bff/Functions.h>

IMPBFF_BEGIN_NAMESPACE

void Functions::shift(double value, std::vector<double> &x)
{
    double ts = -value;
    int tsi = (int) ts;
    double tsf = ts - tsi;

    // shift the values contained in the temporary the vector the integer tsi
    Functions::roll(tsi, x);

    // make a copy and roll it by tsi + 1
    std::vector<double> x_(x);
    Functions::roll(tsi + 1, x_);

    // make the
    // scale the values contained in the temporary y vector by tsf
    for (unsigned int i = 0; i < x.size(); i++) {
        x[i] = x[i] * tsf + x_[i] * (1. - tsf);
    }
}

void Functions::roll(int value, std::vector<double> &y)
{
    if (value > 0) {
        std::rotate(y.begin(), y.begin() + value, y.end());
    } else {
        value = value * (-1);
        std::rotate(y.rbegin(), y.rbegin() + value, y.rend());
    }
}

void Functions::copy_vector_to_array(std::vector<double> &v, double *out, int nout)
{
    for (int i = 0; i < nout; i++) {
        out[i] = v[i];
    }
}

void Functions::copy_array_to_vector(double *in, int nin, std::vector<double> &v)
{
    v.resize(nin);
    for (int i = 0; i < nin; i++) {
        v[i] = in[i];
    }
}

void Functions::copy_vector_to_array(std::vector<double> &v, double **out, int *nout)
{
    *out = (double *) malloc(v.size() * sizeof(double));
    *nout = v.size();
    copy_vector_to_array(v, *out, *nout);
}

void Functions::copy_two_vectors_to_interleaved_array(
        std::vector<double> &v1,
        std::vector<double> &v2,
        double **out, int *nout
)
{
    if (v1.size() == v2.size()) {
        int n = v1.size() + v2.size();
        auto r = (double *) malloc(n * sizeof(double));
        for (unsigned int i = 0; i < v1.size(); i++) {
            r[2 * i + 0] = v1[i];
            r[2 * i + 1] = v2[i];
        }
        *out = r;
        *nout = n;
    }
}

std::vector<double> Functions::diff(std::vector<double> v)
{
    std::vector<double> dx(v.size() - 1);
    for (size_t i = 0; i < dx.size(); i++) {
        dx[i] = (v[i + 1] - v[i]);
    }
    return dx;
}

uint64_t Functions::get_time()
{
    // the birth is the current time stored as a long
    auto now = std::chrono::system_clock::now();
    auto now_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now);
    auto epoch = now_ms.time_since_epoch();
    auto value = std::chrono::duration_cast<std::chrono::milliseconds>(epoch);
    return value.count();
}

bool Functions::bson_iter_skip(bson_iter_t *iter, std::vector<std::string> *skip)
{
    for (auto &sk : *skip) {
        if (strcmp(bson_iter_key(iter), sk.c_str()) == 0) {
            return true;
        }
    }
    return false;
}

void Functions::add_documents(bson_t *src, bson_t *dst, std::vector<std::string> skip)
{
    bson_iter_t iter;
    if (bson_iter_init(&iter, src)) {
        while (bson_iter_next(&iter)) {
            if (!Functions::bson_iter_skip(&iter, &skip)) {
                BSON_APPEND_VALUE(dst, bson_iter_key(&iter), bson_iter_value(&iter));
            }
        }
    }
}

IMPBFF_END_NAMESPACE
