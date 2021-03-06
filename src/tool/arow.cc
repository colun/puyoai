#include <algorithm>
#include <vector>

#include <dirent.h>

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gtest/gtest.h>
#include <SDL_image.h>

#include "base/strings.h"
#include "capture/color.h"
#include "core/real_color.h"
#include "gui/unique_sdl_surface.h"
#include "gui/util.h"
#include "recognition/arow.h"

DECLARE_string(testdata_dir);

using namespace std;

// Before running this trainer, run
// $ tar zxvf recognition.tgz
// in testdata/images directory.

// ----------------------------------------------------------------------

int main()
{
    const int WIDTH = 16;
    const int HEIGHT = 16;

    Arow arows[NUM_RECOGNITION];
    vector<vector<double>> features[NUM_RECOGNITION];

    const pair<string, RecognitionColor> training_testcases[] = {
        make_pair((FLAGS_testdata_dir + "/images/puyo/red.png"), RecognitionColor::RED),
        make_pair((FLAGS_testdata_dir + "/images/puyo/blue.png"), RecognitionColor::BLUE),
        make_pair((FLAGS_testdata_dir + "/images/puyo/yellow.png"), RecognitionColor::YELLOW),
        make_pair((FLAGS_testdata_dir + "/images/puyo/green.png"), RecognitionColor::GREEN),
        make_pair((FLAGS_testdata_dir + "/images/puyo/purple.png"), RecognitionColor::PURPLE),
        make_pair((FLAGS_testdata_dir + "/images/puyo/empty.png"), RecognitionColor::EMPTY),
        make_pair((FLAGS_testdata_dir + "/images/puyo/ojama.png"), RecognitionColor::OJAMA),
        make_pair((FLAGS_testdata_dir + "/images/puyo/zenkeshi.png"), RecognitionColor::ZENKESHI),

        make_pair((FLAGS_testdata_dir + "/images/puyo/red-blur.png"), RecognitionColor::RED),
        make_pair((FLAGS_testdata_dir + "/images/puyo/blue-blur.png"), RecognitionColor::BLUE),
        make_pair((FLAGS_testdata_dir + "/images/puyo/yellow-blur.png"), RecognitionColor::YELLOW),
        make_pair((FLAGS_testdata_dir + "/images/puyo/green-blur.png"), RecognitionColor::GREEN),
        make_pair((FLAGS_testdata_dir + "/images/puyo/purple-blur.png"), RecognitionColor::PURPLE),
        make_pair((FLAGS_testdata_dir + "/images/puyo/empty-blur.png"), RecognitionColor::EMPTY),
        make_pair((FLAGS_testdata_dir + "/images/puyo/ojama-blur.png"), RecognitionColor::OJAMA),
        make_pair((FLAGS_testdata_dir + "/images/puyo/zenkeshi-blur.png"), RecognitionColor::ZENKESHI),

        make_pair((FLAGS_testdata_dir + "/images/puyo/red-actual.png"), RecognitionColor::RED),
        make_pair((FLAGS_testdata_dir + "/images/puyo/blue-actual.png"), RecognitionColor::BLUE),
        make_pair((FLAGS_testdata_dir + "/images/puyo/yellow-actual.png"), RecognitionColor::YELLOW),
        make_pair((FLAGS_testdata_dir + "/images/puyo/green-actual.png"), RecognitionColor::GREEN),
        make_pair((FLAGS_testdata_dir + "/images/puyo/purple-actual.png"), RecognitionColor::PURPLE),
        make_pair((FLAGS_testdata_dir + "/images/puyo/empty-actual.png"), RecognitionColor::EMPTY),
        make_pair((FLAGS_testdata_dir + "/images/puyo/ojama-actual.png"), RecognitionColor::OJAMA),
        make_pair((FLAGS_testdata_dir + "/images/puyo/zenkeshi-actual.png"), RecognitionColor::ZENKESHI),
    };

    // Read training testcases.
    for (const auto& testcase: training_testcases) {
        const string& filename = testcase.first;
        const RecognitionColor color = testcase.second;

        cout << "Opening... " << filename << endl;
        UniqueSDLSurface surf(makeUniqueSDLSurface(IMG_Load(filename.c_str())));
        CHECK(surf.get());

        for (int x = 0; (x + 1) * WIDTH <= surf->w; ++x) {
            for (int y = 0; (y + 1) * HEIGHT <= surf->h; ++y) {
                int pos = 0;
                vector<double> fs(WIDTH * HEIGHT * 3);
                for (int yy = 0; yy < HEIGHT; ++yy) {
                    for (int xx = 0; xx < WIDTH; ++xx) {
                        std::uint32_t c = getpixel(surf.get(), x * WIDTH + xx, y * HEIGHT + yy);
                        std::uint8_t r, g, b;
                        SDL_GetRGB(c, surf->format, &r, &g, &b);
                        fs[pos++] = r;
                        fs[pos++] = g;
                        fs[pos++] = b;
                    }
                }

                CHECK(pos == WIDTH * HEIGHT * 3);
                features[static_cast<int>(color)].push_back(std::move(fs));
            }
        }
    }

    for (const auto& f : features) {
        cout << "SIZE = " << f.size() << endl;
    }

    // training
    for (int times = 0; times < 300; ++times) {
        for (int i = 0; i < NUM_RECOGNITION; ++i) {
            for (int j = 0; j < static_cast<int>(features[i].size()); ++j) {
                if ((j & 0xF) == 0)
                    continue;

                for (int k = 0; k < NUM_RECOGNITION; ++k) {
                    arows[k].update(features[i][j], i == k ? 1 : -1);
                }
            }
        }

        cout << "training " << times << ": done" << endl;
    }

    // test by all
    int num = 0;
    int fail = 0;
    for (int i = 0; i < NUM_RECOGNITION; ++i) {
        for (int j = 0; j < static_cast<int>(features[i].size()); ++j) {
            ++num;
            double vs[NUM_RECOGNITION] {};
            for (int k = 0; k < NUM_RECOGNITION; ++k) {
                vs[k] = arows[k].margin(features[i][j]);
            }

            int result = std::max_element(vs, vs + NUM_RECOGNITION) - vs;
            if (i != result) {
                cout << "fail: " << i << " " << j << " -> " << result << endl;
                for (int k = 0; k < NUM_RECOGNITION; ++k)
                    cout << vs[k] << ' ';
                cout << endl;
                ++fail;
            }
        }
    }

    cout << "num = " << num << endl;
    cout << "fail = " << fail << endl;

    arows[0].save("red.arow");
    arows[1].save("blue.arow");
    arows[2].save("yellow.arow");
    arows[3].save("green.arow");
    arows[4].save("purple.arow");
    arows[5].save("empty.arow");
    arows[6].save("ojama.arow");
    arows[7].save("zenkeshi.arow");

    return 0;
}
