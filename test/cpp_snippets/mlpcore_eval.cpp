// Evaluate a tttrlib.neural_net JSON model with imp.bff's vendored headers only
// (internal/json.h + internal/MlpCore.h) -- the proof that a network trained
// in tttrlib runs, and differentiates, inside bff without linking tttrlib.
//
//   mlpcore_eval <model.json | model.onnx> <X.txt>
// X.txt: first line "n_rows n_cols", then row-major doubles. Prints one line
// per row with the outputs, then a line "fd_check <max|dL/dparams - FD|>" for
// the loss L = sum(y) over the batch, then "dx_check <max|dL/dx - FD|>".
#define TTTRLIB_MLPCORE_NAMESPACE IMP::bff::internal
#include "IMP/bff/internal/json.h"
#include "IMP/bff/internal/MlpCore.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>

namespace mc = IMP::bff::internal::mlpcore;
using IMP::bff::internal::MlpModel;

int main(int argc, char** argv) {
    if (argc < 3) return 2;
    // .json: the tttrlib.neural_net document; .onnx: straight from any exporter
    const std::string path(argv[1]);
    MlpModel m;
    if (path.size() > 5 && path.compare(path.size() - 5, 5, ".onnx") == 0) {
        std::ifstream f(path, std::ios::binary);
        std::stringstream ss;
        ss << f.rdbuf();
        m = mc::model_from_onnx(ss.str());
    } else {
        std::ifstream fj(path);
        nlohmann::json j;
        fj >> j;
        m = mc::model_from_json(j);
    }

    std::ifstream fx(argv[2]);
    int n_rows, n_cols;
    fx >> n_rows >> n_cols;
    std::vector<double> X(static_cast<size_t>(n_rows) * n_cols);
    for (auto& v : X) fx >> v;

    std::vector<double> y, d1, d2;
    mc::model_predict(m, X.data(), n_rows, 0, nullptr, y, d1, d2);
    const int n_out = m.n_outputs();
    for (int r = 0; r < n_rows; ++r) {
        for (int k = 0; k < n_out; ++k) std::printf("%.17g ", y[static_cast<size_t>(r) * n_out + k]);
        std::printf("\n");
    }

    // L = sum(y): dL/dy = 1
    std::vector<double> ones(y.size(), 1.0), dp, dX, dV;
    mc::model_backward(m, X.data(), n_rows, nullptr, ones.data(), nullptr, nullptr, dp, dX, dV);
    auto loss = [&](const MlpModel& mm, const std::vector<double>& x) {
        std::vector<double> a, b, c;
        mc::model_predict(mm, x.data(), n_rows, 0, nullptr, a, b, c);
        double s = 0; for (double v : a) s += v; return s;
    };
    const double h = 1e-5;
    std::vector<double> p; mc::flatten(m.layers, p);
    double worst = 0;
    for (size_t i = 0; i < p.size(); ++i) {
        MlpModel mp = m, mm = m; auto pp = p, pm = p; pp[i] += h; pm[i] -= h;
        mc::unflatten(mp.layers, pp.data(), pp.size()); mc::unflatten(mm.layers, pm.data(), pm.size());
        const double fd = (loss(mp, X) - loss(mm, X)) / (2 * h);
        worst = std::max(worst, std::abs(fd - dp[i]));
    }
    std::printf("fd_check %.3g\n", worst);
    double worst_x = 0;
    for (size_t i = 0; i < X.size(); ++i) {
        auto xp = X, xm = X; xp[i] += h; xm[i] -= h;
        worst_x = std::max(worst_x, std::abs((loss(m, xp) - loss(m, xm)) / (2 * h) - dX[i]));
    }
    std::printf("dx_check %.3g\n", worst_x);
    return 0;
}
