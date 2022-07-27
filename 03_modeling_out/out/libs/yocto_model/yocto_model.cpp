//
// Implementation for Yocto/Model
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_model.h"

#include <yocto/yocto_sampling.h>

#include "ext/perlin-noise/noise1234.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::string;
using std::vector;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR EXAMPLE OF PROCEDURAL MODELING
// -----------------------------------------------------------------------------
namespace yocto {

float noise(const vec3f& p) { return ::noise3(p.x, p.y, p.z); }
vec2f noise2(const vec3f& p) {
  return {noise(p + vec3f{0, 0, 0}), noise(p + vec3f{3, 7, 11})};
}
vec3f noise3(const vec3f& p) {
  return {noise(p + vec3f{0, 0, 0}), noise(p + vec3f{3, 7, 11}),
      noise(p + vec3f{13, 17, 19})};
}
float fbm(const vec3f& p, int octaves) {
  auto sum    = 0.0f;
  auto weight = 1.0f;
  auto scale  = 1.0f;
  for (auto octave = 0; octave < octaves; octave++) {
    sum += weight * fabs(noise(p * scale));
    weight /= 2;
    scale *= 2;
  }
  return sum;
}
float turbulence(const vec3f& p, int octaves) {
  auto sum    = 0.0f;
  auto weight = 1.0f;
  auto scale  = 1.0f;
  for (auto octave = 0; octave < octaves; octave++) {
    sum += weight * fabs(noise(p * scale));
    weight /= 2;
    scale *= 2;
  }
  return sum;
}
float ridge(const vec3f& p, int octaves) {
  auto sum    = 0.0f;
  auto weight = 0.5f;
  auto scale  = 1.0f;
  for (auto octave = 0; octave < octaves; octave++) {
    sum += weight * (1 - fabs(noise(p * scale))) * (1 - fabs(noise(p * scale)));
    weight /= 2;
    scale *= 2;
  }
  return sum;
}

void add_polyline(shape_data& shape, const vector<vec3f>& positions,
    const vector<vec4f>& colors, float thickness = 0.0001f) {
  auto offset = (int)shape.positions.size();
  shape.positions.insert(
      shape.positions.end(), positions.begin(), positions.end());
  shape.colors.insert(shape.colors.end(), colors.begin(), colors.end());
  shape.radius.insert(shape.radius.end(), positions.size(), thickness);
  for (auto idx = 0; idx < positions.size() - 1; idx++) {
    shape.lines.push_back({offset + idx, offset + idx + 1});
  }
}

void sample_shape(vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, const shape_data& shape, int num) {
  auto triangles  = shape.triangles;
  auto qtriangles = quads_to_triangles(shape.quads);
  triangles.insert(triangles.end(), qtriangles.begin(), qtriangles.end());
  auto cdf = sample_triangles_cdf(triangles, shape.positions);
  auto rng = make_rng(19873991);
  for (auto idx = 0; idx < num; idx++) {
    auto [elem, uv] = sample_triangles(cdf, rand1f(rng), rand2f(rng));
    auto q          = triangles[elem];
    positions.push_back(interpolate_triangle(
        shape.positions[q.x], shape.positions[q.y], shape.positions[q.z], uv));
    normals.push_back(normalize(interpolate_triangle(
        shape.normals[q.x], shape.normals[q.y], shape.normals[q.z], uv)));
    if (!texcoords.empty()) {
      texcoords.push_back(interpolate_triangle(shape.texcoords[q.x],
          shape.texcoords[q.y], shape.texcoords[q.z], uv));
    } else {
      texcoords.push_back(uv);
    }
  }
}

void make_terrain(shape_data& shape, const terrain_params& params) {
  // YOUR CODE GOES HERE
}

void make_displacement(shape_data& shape, const displacement_params& params) {
  // YOUR CODE GOES HERE
}

void make_hair(
    shape_data& hair, const shape_data& shape, const hair_params& params) {
  // YOUR CODE GOES HERE
}

void make_grass(scene_data& scene, const instance_data& object,
    const vector<instance_data>& grasses, const grass_params& params) {
  // YOUR CODE GOES HERE
}

}  // namespace yocto
