// adapted https://stackoverflow.com/questions/10990618/calculate-rolling-moving-average-in-c
template <typename T, typename Total>
class MovingAverage {

public:
    MovingAverage(size_t window_size){
        samples_.resize(window_size);
        this->window_size = window_size;
    }
    MovingAverage& operator()(T sample)
    {
        total_ += sample;
        if (num_samples_ < window_size)
            samples_[num_samples_++] = sample;
        else
        {
            T& oldest = samples_[num_samples_++ % window_size];
            total_ -= oldest;
            oldest = sample;
        }
        return *this;
    }
    operator double() const { return total_ / std::min(num_samples_, window_size); }
private:
    size_t window_size;
    std::vector<T> samples_;
    size_t num_samples_{0};
    Total total_{0};
};