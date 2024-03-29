import typing
import numpy as np

import IMP.bff

from typing import NamedTuple

MAX_FLOAT = np.finfo(np.float64).max
MIN_FLOAT = np.finfo(np.float64).min


class DataSettings(NamedTuple):
    data: np.ndarray = None
    data_noise: np.ndarray = None
    time_axis: np.array = None
    irf_histogram: np.array = None
    acquisition_time: float = MAX_FLOAT


class LifetimeSettings(NamedTuple):
    lifetime_spectrum: np.array = np.array([], dtype=np.float64)
    use_amplitude_threshold: bool = False
    abs_lifetime_spectrum: bool = True
    amplitude_threshold: float = MIN_FLOAT


class ConvolutionSettings(NamedTuple):
    convolution_method: int = IMP.bff.DecayConvolution.FAST
    excitation_period: float = 100.0
    irf_shift_channels: float = 0.0
    irf_background_counts: float = 0.0
    start: int = 0
    stop: int = -1
    active: bool = True


class PatternSettings(NamedTuple):
    constant_offset: float = 0.0
    pattern: IMP.bff.DecayCurve = None
    pattern_fraction: float = 0.0
    start: int = 0
    stop: int = -1
    active: bool = True


class PileupSettings(NamedTuple):
    pile_up_model: str = "coates"
    repetition_rate: float = 10.0
    instrument_dead_time: float = 120.0
    start: int = 0
    stop: int = -1
    active: bool = True


class ScalingSettings(NamedTuple):
    constant_offset: float = 0.0
    start: int = 0
    stop: int = -1
    active: bool = True


class LinearizationSettings(NamedTuple):
    linearization_table: IMP.bff.DecayCurve = None
    start: int = 0
    stop: int = -1
    active: bool = True


class ScoreSettings(NamedTuple):
    score_type: str = "default"
    start: int = 0
    stop: int = -1


class Decay(object):

    # These attributes are skipped when saving to JSON
    _ignore_attributes = [
        'irf_histogram',
        # 'tttr_data',
        # 'tttr_irf',
        'linearization_table',
        'corrected_irf',
        'x',
        'y',
        'data_weights',
        'number_of_photons'
    ]

    _attributes = {
        #VERSION
        #####################
        # TODO
        # 'version': ("", "get_version"),
        # Data
        ############################
        'data': ("_data", "y"),
        'data_noise': ("_data", "ey"),
        'model': ("_model", "y"),
        'y': ("_model", "y"),
        'x': ("_model", "x"),
        'acquisition_time': ("_data", "acquisition_time"),
        'time_axis': ("_data", "x"),
        'irf_histogram': ("_irf", "y"),
        'irf': ("_irf", "y"),
        # 'tttr_data': ("_data", None, "set_tttr"),
        # 'tttr_irf': ("_irf", None, "set_tttr"),
        # LIFETIME SPECTRUM
        ############################
        'lifetime_spectrum': ("lifetime_handler", "get_lifetime_spectrum", "set_lifetime_spectrum"),
        'use_amplitude_threshold': ("lifetime_handler", "use_amplitude_threshold"),
        'abs_lifetime_spectrum': ("lifetime_handler", "abs_lifetime_spectrum"),
        'amplitude_threshold': ("lifetime_handler", "amplitude_threshold"),
        # CONVOLUTION
        ############################
        'corrected_irf': ("decay_convolution.corrected_irf", "y"),
        'convolution_range': ("decay_convolution", "get_range", "set_range"),
        'convolution_method': ("decay_convolution", "convolution_method"),
        'excitation_period': ("decay_convolution", "excitation_period"),
        'irf_shift_channels': ("decay_convolution", "irf_shift_channels"),
        'irf_background_counts': ("decay_convolution", "irf_background_counts"),
        # Scatter
        ############################
        'scatter_fraction': ("decay_scatter", "pattern_fraction"),
        # BACKGROUND
        ############################
        'constant_offset': ("", "constant_background"),
        'pattern_fraction': ("decay_background", "pattern_fraction"),
        # PILE-UP
        #############################
        'pile_up_model': ("decay_pileup", "pile_up_model"),
        'instrument_dead_time': ("decay_pileup", "instrument_dead_time"),
        'use_pile_up_correction': ("decay_pileup", "active"),
        # SCALE
        ############################
        'scale_model_to_data': ("decay_scale", "active"),
        'number_of_photons': ("decay_scale", "number_of_photons"),
        # LINEARIZATION
        ############################
        'linearization_table': ("decay_linearization.data", "y"),
        'linearization': ("decay_linearization.data", "y"),
        'use_linearization': ("decay_linearization", "active"),
        # SCORE
        ###########################
        'score_type': ("decay_score", "score_type"),
        'score_range': ('', 'get_score_range', 'set_score_range'),
        'score': ("decay_score", "score"),
        'weighted_residuals': ("decay_score", "weighted_residuals")
    }

    _methods = {
        'get_mean_lifetime': ("decay_convolution", "get_mean_lifetime"),
        'get_score': ("decay_score", "get_score")
    }

    @staticmethod
    def make_decay_curves(
            data: np.ndarray = None,
            data_noise: np.ndarray = None,
            time_axis: np.array = None,
            irf_histogram: np.array = None,
            # tttr_data: tttrlib.TTTR = None,
            # tttr_irf: tttrlib.TTTR = None,
            # tttr_micro_time_coarsening: int = 1,
            acquisition_time: float = MAX_FLOAT,
    ) -> typing.Tuple[IMP.bff.DecayCurve, IMP.bff.DecayCurve]:
        # Data
        if time_axis is None:
            time_axis = np.array([], np.float64)
        if data is None:
            data = np.array([], np.float64)
        # IRF
        if irf_histogram is None:
            irf_histogram = np.zeros_like(data)
            if len(irf_histogram) > 0:
                irf_histogram[0] = 1.0
        decay = IMP.bff.DecayCurve(x=time_axis, y=data, acquisition_time=acquisition_time)
        irf = IMP.bff.DecayCurve(x=time_axis, y=irf_histogram)
        if data_noise is not None:
            decay.set_ey(data_noise)
        return decay, irf

    def __init__(
        self,
        data_settings=None,  # type: DataSettings
        lifetime_settings=None,  # type: LifetimeSettings
        convolution_settings=None,  # type: ConvolutionSettings
        scatter_settings=None,  # type: PatternSettings
        background_settings=None,  # type: PatternSettings
        pileup_settings=None,  # type: PileupSettings
        scaling_settings=None,  # type: ScalingSettings
        linearization_settings=None,  # type: LinearizationSettings
        score_settings=None  # type: ScoreSettings
    ):
        super().__init__()
        if data_settings is None:
            data_settings = DataSettings()
        if lifetime_settings is None:
            lifetime_settings = LifetimeSettings()
        if convolution_settings is None:
            convolution_settings = ConvolutionSettings()
        if background_settings is None:
            background_settings = PatternSettings()
        if scatter_settings is None:
            scatter_settings = PatternSettings()
        if pileup_settings is None:
            pileup_settings = PileupSettings()
        if linearization_settings is None:
            linearization_settings = LinearizationSettings()
        if scaling_settings is None:
            scaling_settings = ScalingSettings()
        if score_settings is None:
            score_settings = ScoreSettings()

        # init decay curves
        data, irf = self.make_decay_curves(*data_settings)
        model = IMP.bff.DecayCurve(data.x)

        lifetime_handler = IMP.bff.DecayLifetimeHandler(*lifetime_settings)
        decay_convolution = IMP.bff.DecayConvolution(lifetime_handler, irf, *convolution_settings)
        decay_scatter = IMP.bff.DecayPattern(*scatter_settings)
        decay_scatter.data = irf
        decay_background = IMP.bff.DecayPattern(*background_settings)
        decay_pileup = IMP.bff.DecayPileup(data, *pileup_settings)
        decay_linearization = IMP.bff.DecayLinearization(*linearization_settings)
        decay_scale = IMP.bff.DecayScale(data, *scaling_settings)
        decay_score = IMP.bff.DecayScore(model, data, *score_settings)

        self._data, self._irf, self._model = data, irf, model
        self.lifetime_handler = lifetime_handler
        self.decay_convolution = decay_convolution
        self.decay_scatter = decay_scatter
        self.decay_pileup = decay_pileup
        self.decay_linearization = decay_linearization
        self.decay_scale = decay_scale
        self.decay_background = decay_background
        self.decay_score = decay_score

        self._decay_modifier = [
            self.decay_convolution,
            self.decay_scatter,
            # self.decay_pileup, # TODO: seems cn
            self.decay_linearization,
            self.decay_scale,
            self.decay_background,
        ]

    def load(self, *args, **kwargs):
        super().load(*args, **kwargs)
        self._model.x = self._data.x
        self._irf.x = self._data.x

    def update(self):
        dc = self
        for dm in self._decay_modifier:
            dm.add(dc._model)

    @property
    def score(self, *args, **kwargs):
        self.update()
        return self.decay_score.get_score(*args, **kwargs)

    @property
    def mean_lifetime(self):
        return self.decay_convolution.get_mean_lifetime(self._data)

    @property
    def constant_background(self):
        return self.decay_scale.constant_background

    @constant_background.setter
    def constant_background(self, v):
        self.decay_scale.constant_background = v
        self.decay_background.constant_offset = v

    def get_score_range(self):
        return self.decay_score.get_range(self._data)

    def set_score_range(self, start_stop):
        for dm in self._decay_modifier:
            dm.set_range(start_stop)
        self.decay_score.set_range(start_stop)

