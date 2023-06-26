#include "DatasetOperations.h"

//using std::move;
//using std::size_t;
//using std::string;
//using std::vector;

constexpr size_t NORMALIZE_DENOM = 255;

template <template <typename...> class Container, typename Image>
void normalize(Container<Image>& images)
{
    size_t image_count = images.size();
    size_t pixel_count_per_image = images[0].size();
    for (size_t i = 0; i < image_count; ++i)
    {
        for (size_t j = 0; j < pixel_count_per_image; ++j)
        {
            images[i][j] /= NORMALIZE_DENOM;
        }
    }
}

vector<vector<float>> loadMnistTestImages(const std::string& dataset_dir, const size_t& test_limit)
{
    mnist::MNIST_dataset<vector, vector<float>, uint8_t> dataset = mnist::read_dataset<vector, vector, float, uint8_t>(dataset_dir, 1, test_limit);
    normalize(dataset.test_images);

    return dataset.test_images;
}

vector<unsigned char> loadMnistTestLabels(const std::string& dataset_dir, const size_t& test_limit)
{
    mnist::MNIST_dataset<vector, vector<float>, uint8_t> dataset = mnist::read_dataset<vector, vector, float, uint8_t>(dataset_dir, 1, test_limit);

    return dataset.test_labels;
}
