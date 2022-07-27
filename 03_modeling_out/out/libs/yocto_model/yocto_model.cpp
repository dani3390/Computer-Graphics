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

#include "iostream"
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


vec2f hash2(vec2f p) {
  vec2f q = vec2f(
      {dot(p, vec2f({127.1f, 311.7f})), dot(p, vec2f({269.5f, 183.3f}))});
  vec2f x = vec2f({sin(q.x) * 43758.5453f, sin(q.y) * 43758.5453f});
  return x - vec2f({floor(x.x), floor(x.y)});
  // voronoise position.x e y
}

vec3f hash3(vec3f p) {
  vec3f q = vec3f({dot(p, vec3f({127.1f, 311.7f, 1.0f})),
      dot(p, vec3f({269.5f, 183.3f, 1.0f})),
      dot(p, vec3f({419.2f, 371.9f, 1.0f}))});
  vec3f x = vec3f({sin(q.x) * 43758.5453f, sin(q.y) * 43758.5453f, sin(q.z) * 43758.5453f});
  return x - vec3f({floor(x.x), floor(x.y), floor(x.z)});
  // voronoise position.x e y
}

float func2(vec2f position) {
  vec2f p = vec2f({floor(position.x), floor(position.y)});
  vec2f f = position - p;

  vec2f mb;
  vec2f mr;
  auto  rng = make_rng(100);

  float res = 8.0f;
  for (int j = -1; j <= 1; j++)
    for (int i = -1; i <= 1; i++) {
      vec2f b = vec2f({float(i), float(j)});

      vec2f o = hash2(p + b);

      vec2f r = b + rand2f(rng) - f;  // rand2f(rng)
      float d = dot(r, r);

      if (d < res) {
        res = d;
        mr  = r;
        mb  = b;
      }
    }

  res = 8.0f;
  for (int j = -2; j <= 2; j++)
    for (int i = -2; i <= 2; i++) {
      vec2f b = mb + vec2f{float(i), float(j)};
      vec2f o = hash2(p + b);
      vec2f r = b + rand2f(rng) - f;
      //  std::cout << r.x << " " << r.y << std::endl;
      float d = dot(0.5 * (mr + r), normalize(r - mr));

      res = min(res, d);
    }

  return res;
}

vec3f cell(vec2f position) {
  vec2f p = vec2f{floor(position.x), floor(position.y)};
  vec2f f = position - p;

  vec2f mb;
  vec2f mr;
  auto  rng = make_rng(100);

  float res = 8.0;
  for (int j = -1; j <= 1; j++)
    for (int i = -1; i <= 1; i++) {
      vec2f b = vec2f({float(i), float(j)});

      vec2f o = hash2(p + b);

      vec2f r = b + o - f;  // rand2f(rng)
      float d = dot(r, r);

      if (d < res) {
        res = d;
        mr  = r;
        mb  = b;
      }
    }

  res = 8.0f;
  for (int j = -2; j <= 2; j++)
    for (int i = -2; i <= 2; i++) {
      vec2f b = mb + vec2f{float(i), float(j)};
      vec2f o = hash2(p + b);
      vec2f r = b + o - f;
      //  std::cout << r.x << " " << r.y << std::endl;
      if (dot(mr - r, mr - r) > 0.00001f)
        res = min(res, dot(0.5 * (mr + r), normalize(r - mr)));
    }

  return vec3f{(res, mr.x, mr.y)};
}

vec3f cell2(vec2f p) { //, float u, float v
  vec3f d = cell(p);

  vec3f col = d.x * (0.5f + 0.5f * sin(64.0f * d.x)) *
              vec3f({1.0f, 1.0f, 1.0f});

  col = interpolate_line(
      vec3f({1.0f, 0.6f, 0.0f}), col, smoothstep(0.04f, 0.07f, d.x));

 float dd = length(vec2f({d.y, d.z}));

  col = interpolate_line(
      vec3f({1.0f, 0.6f, 0.1f}), col, smoothstep(0.0f, 0.12f, dd));
 
  col += vec3f({1.0f, 0.6f, 0.1f}) * (1.0f - smoothstep(0.0f, 0.04f, dd));

  return col;
}


float extra2(vec2f p, float u, float v) 
{
  float d = func2(p);

  return 1.0f - smoothstep(u, v, d);
}

float voronoise(vec3f position, float u, float v) {
  float k = 1.0 + 63.0 * pow(1.0 - v, 6.0);

  
  vec3f i = {floor(position.x), floor(position.y), floor(position.z)};
  vec3f f = position - i;
  vec3f a = vec3f({0.0f, 0.0f, 0.0f});
  
  float va = 0.0f;
  float wt = 0.0f;

  for (int z = -2; z <= 2; z++)
  for (int y = -2; y <= 2; y++)
    for (int x = -2; x <= 2; x++) {
      vec3f g  = {float(x), float(y), float(z)};
      vec3f o  = hash3(i + g) * vec3f({u, u, 1.0f});
      vec3f d  = g - f + o.x * o.y * o.z;
      float da = dot(d, d);
      float w  = pow(1.0f - smoothstep(0.0f, 1.414f, sqrt(da)), k);
      va += o.z * w;
      wt += w;
    }

  float color = va / wt;

  return color;
  //      vec3f next = positions[k] + t * normal +
  // noise3(positions[k] * params.scale) * params.strength;
}

void make_terrain(shape_data& shape, const terrain_params& params) {
  for (int i = 0; i < shape.positions.size(); i++) {
    vec3f position = shape.positions[i];
    position += shape.normals[i] *
                ridge(position * params.scale, params.octaves) * params.height *
                (1.f - length(position - params.center) / params.size);
    float perc = position.y / params.height;
    if (perc <= 0.3f)
      shape.colors.push_back(params.bottom);
    else if (perc > 0.6f)
      shape.colors.push_back(params.top);
    else
      shape.colors.push_back(params.middle);
    shape.positions[i] = position;
  }
  compute_normals(shape.normals, shape);
}

float smoothvoronoi(vec2f position, float u) {


      vec2f l = {floor(position.x), floor(position.y)};
  vec2f f = position - l;

  float res = 0.0;
  for (int j = -1; j <= 1; j++)
    for (int i = -1; i <= 1; i++) {
      vec2f b = vec2f({float(i), float(j)});
      vec2f r = b - f + hash2(l + b);
      float d = dot(r, r);

      res += 1.0f / pow(d, 8.0f);
    }
  return pow(u / res, u / 16.0f);


}

/* vec2f eval_cell(vec2f gij, vec2f at) { 
  vec2f res  = vec2f({0.0f, 0.0f});
  int  gidx = coord2index(gij, grid_size);
  for (int i = 0; i < nb_Kernel; i++) {
    // compute the kernel coordinate id
    int kid = gidx * nb_Kernel * 2 + i * 2;
    // compute the kernel coordinate in the buffer
    ivec2 c1 = index2coord(kid, resx), c2 = index2coord(kid + 1, resx);
    // fetch the kernel k
    vec4 k1 = texelFetch(iChannel0, c1, 0), k2 = texelFetch(iChannel0, c2, 0);
    // flag for the interface
    float phi = optimPhase ? k1.w : 0., filtering = filterPhase ? 1. : 0.;
    // computing the impacte of the phasor kernel
    res += phasorkernel(at - k1.xy, k2.w, k2.xy, phi, filtering);
  }
  return res;
}


vec2f phasor(vec2f gij, vec2f ij) {
  vec2f res = vec2f({0.0f, 0.0f});
  // loop on the neighbor cells
  for (int i = -2; i <= 2; i++)
    for (int j = -2; j <= 2; j++)
      // evaluation of the cell
      res += eval_cell(
          gij + vec2f({float(i), float(j)}), ij - vec2f({float(i), float(j)}));

  return res;
}*/

void make_displacement(shape_data& shape, const displacement_params& params) {

    float cellsize = 1.0f / 30.0f; 
  for (int i = 0; i < shape.positions.size(); i++) {
    
    vec3f position = shape.positions[i];
    vec3f last_p   = vec3f(position); 
   
    vec3f pos = position / params.height / cellsize;
    vec2f ij  = vec2f({pos.x, pos.y}) - vec2f({floor(pos.x), floor(pos.y)});
    vec2f gij = vec2f({floor(pos.x + 1.0f), floor(pos.y + 1.0f)});
   

    if (params.extra2) {

        position += shape.normals[i] *
                  extra2(vec2f({position.x, position.y}) * params.scale,
                      params.u, params.v) *
                  params.height;

    } 
    else if (params.voronoise) {
      position += shape.normals[i] *
          voronoise(vec3f({position.x, position.y, position.z}) * params.scale,
              params.u, 1.0f) *
                  params.height;
    
    } 
    else if (params.smooth) {
      position +=
          shape.normals[i] *
                  smoothvoronoi(vec2f({position.x, position.y}) * params.scale,
                      params.u) *
          params.height;
    
    } 
    else if (params.cellnoise) {

      position +=
          shape.normals[i] *
          voronoise(vec3f({position.x, position.y, position.z}) * params.scale,
              params.u, params.v) *
          params.height;
    

    }
    
    else {
      position += shape.normals[i] *
                  (turbulence(position * params.scale, params.octaves) *
                      params.height);
    
    }

    shape.colors.push_back(interpolate_line(params.bottom, params.top, distance(position, last_p) / params.height));
    shape.positions[i] = position;
  }
  compute_normals(shape.normals, shape);
}

void make_hair(
    shape_data& hair, const shape_data& shape, const hair_params& params) {
  auto shape2 = shape;

  int size = shape2.positions.size();
  sample_shape(
      shape2.positions, shape2.normals, shape2.texcoords, shape2, params.num);
  float t = params.lenght / params.steps;


    auto rng = make_rng(8238328);

  for (int i = size; i < shape2.positions.size(); i++) {
    vector<vec3f> positions;
    vector<vec4f> colors;


      auto val = rand1f(rng);

      if (val < params.density) {
        continue;
      }



    positions.push_back(shape2.positions[i]);
    colors.push_back(params.bottom);


    vec3f normal = shape2.normals[i];
    for (int k = 0; k < params.steps; k++) {

      vec3f next = positions[k] + t * normal +
                   noise3(positions[k] * params.scale) * params.strength;

      next.y -= params.gravity;

      normal = normalize(next - positions[k]); 
      positions.push_back(next);

      colors.push_back(interpolate_line(params.bottom, params.top,
          (distance(next, positions[0]) / params.lenght)));


    }

    colors[params.steps] = params.top;

    add_polyline(hair, positions, colors);
  }

  auto tang = lines_tangents(hair.lines, hair.positions);
  for (int i = 0; i < tang.size(); i++)
    hair.tangents.push_back(vec4f{(tang[i], 0.f)});

 
}

void make_grass(scene_data& scene, const instance_data& object,
    const vector<instance_data>& grasses, const grass_params& params) {
  auto rng = make_rng(198767);

  auto shape = scene.shapes[object.shape];


  sample_shape(
      shape.positions, shape.normals, shape.texcoords, shape, params.num);

  for (int i = 0; i < shape.positions.size(); i++) {
    vec3f         position = shape.positions[i];
    instance_data new_el   = object;

    auto grass = grasses[rand1i(rng, grasses.size())];


    new_el.shape    = grass.shape;
    new_el.material = grass.material;


    new_el.frame.y = shape.normals[i];
    new_el.frame.x = normalize(
        vec3f({1, 0, 0}) -
        dot(vec3f({1, 0, 0}), new_el.frame.y) * new_el.frame.y);
    new_el.frame.z = cross(new_el.frame.x, new_el.frame.y);
    new_el.frame.o = position;

    float rand = 0.9f + rand1f(rng) * 0.1f;

    new_el.frame *= scaling_frame(vec3f({rand, rand, rand}));

    rand = rand1f(rng) * 2 * pif;

    new_el.frame *= rotation_frame(new_el.frame.y, rand);

    rand = 0.1f + rand1f(rng) * 0.1f;

    new_el.frame *= rotation_frame(new_el.frame.z, rand);

    scene.instances.push_back(new_el);
  }
}

}  // namespace yocto
