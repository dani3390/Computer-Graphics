//
// Implementation for Yocto/Grade.
//

//
// LICENSE:
//
// Copyright (c) 2020 -- 2020 Fabio Pellacini
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

#include "yocto_colorgrade.h"

#include <yocto/yocto_cli.h>

#include <yocto/yocto_color.h>
#include <yocto/yocto_sampling.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_color.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/ext/stb_image.h>
#include <yocto/ext/stb_image_write.h>
#include <yocto/ext/tinyexr.h>

#include "stdint.h"


#include "iostream"


using namespace std;


#ifdef GL_ES
precision mediump float;
#endif

float normpdf(float x, float sigma) {
  return 0.39894 * exp(-0.5 * x * x / (sigma * sigma)) / sigma;
}




// -----------------------------------------------------------------------------
// COLOR GRADING FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {



color_image grade_image(const color_image& image, const grade_params& params) {
  auto        img       = image;

  auto img3 = image;
  int  r    = 0;

  auto rng = make_rng(1254664589, 1Ui64);

  int acc = 0, azz = 0;

  auto filename = "apps/imagina.hdr";

  auto error = string{};
  auto img2  = image_data{};
  if (!load_image(filename, img2, error)) {
    cout << error << endl;
  }

  for (auto j = 0; j < img.height; j++) {
    for (auto i = 0; i < img.width; i++) {
      vec2f lalalu = {(i, j)};
      vec2f sizeu   = {img.width, img.height};
      auto pixel = img[{i, j}];
      auto c     = xyz(pixel);

      if (params.exposure != 0) {
        c = c * pow(2, params.exposure);
        // cout << "prova" << endl;
      }

      if (params.filmic == true) {
        c *= 0.6;
        c = (pow(c, 2) * 2.51 + c * 0.03) /
            (pow(c, 2) * 2.43 + c * 0.59 + 0.14);
      }

      if (params.srgb == true) c = pow(c, 1 / 2.2);

      c = clamp(c, 0, 1);

      c = {c.x * params.tint.x, c.y * params.tint.y, c.z * params.tint.z};

      //  if (params.saturation != 0) {
      auto g = (c.x + c.y + c.z) / 3;
      c      = g + (c - g) * (params.saturation * 2);

      // }

      if (params.inverted == true) {
        c = 1.0 - c;
      }

  

      if (params.contrast != 0) {
        c = gain(c, 1 - params.contrast);
      }

      if (params.vignette != 0) {
        auto  vr   = 1 - params.vignette;
        vec2f size = {img.width, img.height};
        vec2f ij   = {i, j};
        auto  r    = length(ij - size / 2) / length(size / 2);
        c          = c * (1 - smoothstep(vr, 2 * vr, r));
      }

      

      if (params.grain != 0) {
        c = c + (rand1f(rng) - 0.5) * params.grain;
      }


      set_pixel(img, i, j, {c.x, c.y, c.z, pixel.w});

      if ((params.vignetterectbt > 0) && (params.vignetterectlr > 0)) {
        vec2f size = {img.width, img.height};
        vec2f ij   = {i, j};

        auto rt = vec2f{0.0, params.vignetterectbt};
        auto rb = vec2f{0.0, params.vignetterectbt};
        auto rl = vec2f{0.0, params.vignetterectlr};
        auto rr = vec2f{0.0, params.vignetterectlr};

        vec2f vu = ij / size;

        auto vt = smoothstep(1 - rt.x, 1 - rt.y, vu.y);

        auto vb = smoothstep(rb.x, rb.y, vu.y);
        auto vl = smoothstep(rl.x, rl.y, vu.x);
        auto vr = smoothstep(1 - rr.x, 1 - rr.y, vu.x);

        c *= (vt * vb * vl * vr);

        set_pixel(img, i, j, {c.x, c.y, c.z, 1.0});
      }

      if (params.gray == true) {
        vec2f size = {img.width, img.height};
        vec2f ij   = {i, j};

        vec2f uv = ij / size;

        auto  pixel = img[{i, j}];
        vec3f m     = xyz(pixel);

        float pruv = 0.333f;
        vec3f x =
            m *
            (pruv *
                2);  // m.x * (pruv * 2) + m.y * (pruv * 2) + m.z * (pruv * 2);

        auto pixel2 = img2[{i, j}];
        auto l      = xyz(pixel2);

        vec3f lolu = l * (pruv);  // l.x * pruv + l.y * pruv + l.z * pruv;

        float avg = dot(x, lolu);

        float vig = 1.0f - distance(uv, vec2f({0.5, 0.5}));

        vec3f bada = {1.0f, 1.0f, 1.0f};

        vec3f lullo = bada * avg * 1.4f * vig;

        set_pixel(img, i, j, {lullo.x, lullo.y, lullo.z, pixel.w});

        vec2f vu = 1 - ij / size;
      }

      if (params.sigma != 0) {
        vec2f size = {img.width, img.height};

               const int mSize = 19;
               const int kSize = (mSize - 1) / 2;
               float     kernel[mSize];
               vec3f     final_colour = vec3f{(0.0f, 0.0f, 0.0f)};

               float sigma = params.sigma;
               float Z     = 0.0;
               for (int q = 0; q <= kSize; ++q) {
                 kernel[kSize + q] = kernel[kSize - q] = normpdf(float(q), sigma);
               }

                  for (int q = 0; q < mSize; ++q) { 
                      Z += kernel[q];
               }


               for (int w = -kSize; w <= kSize; ++w) {
                 for (int m = -kSize; m <= kSize; ++m) {
                   if (j + m <= 1024 - mSize && i + w <= 2048 - mSize) {



                       auto trullu = get_pixel(img, i + w + kSize, j + m + kSize);

                     auto trullulo = xyz(trullu);

                     final_colour += kernel[kSize + m] * kernel[kSize + w] *
                                     trullulo;

                   }

                 }
               }

               auto itru    = final_colour / (Z * Z);

               set_pixel(img, i, j, {itru.x, itru.y, itru.z, 1.0});
      }
    }
  }

  auto red = vec3f{0.393, 0.769, 0.189}, green = vec3f{0.349, 0.686, 0.168},
       blue = vec3f{0.272, 0.534, 0.131};


  mat3f prova({red, green, blue});


         



    if (params.vintage == true) {
       for (auto j = 0; j < img.height; j++) {
           for (auto i = 0; i < img.width; i++) {
                auto pixel = img[{i, j}];
                auto c     = xyz(pixel);
                auto d     = xyz_to_rgb(c);

                auto output = d * prova;

                set_pixel(img, i, j, {output.x, output.y, output.z, pixel.w});
           
                



      }
    }
    
  }


    



     if (params.mosaic != 0) {
         for (auto j = 0; j < img.height; j++) {
            for (auto i = 0; i < img.width; i++) {
             img[{i, j}] = img[{i - i % params.mosaic, j - j % params.mosaic}];


           /* auto pixel = img[{i, j}];
                auto c     = xyz(pixel);
                
                c[i, j] = c[i - i % params.mosaic, j - j % params.mosaic];

                set_pixel(img, i, j, {c.x, c.y, c.z, pixel.w});*/
            }
         
         }


     
     }

     if (params.grid != 0) {
       for (auto j = 0; j < img.height; j++) {
         for (auto i = 0; i < img.width; i++) {
           img[{i, j}] = (0 == i % params.grid || 0 == j % params.grid)
                             ? 0.5f * img[{i, j}]
                             : img[{i, j}];

         

           /* auto pixel = img[{i, j}];
                auto c     = xyz(pixel);

                c[i, j] = c[i - i % params.mosaic, j - j % params.mosaic];

                set_pixel(img, i, j, {c.x, c.y, c.z, pixel.w});*/
         }
       }
     }





  return img;

}

}
// namespace yocto