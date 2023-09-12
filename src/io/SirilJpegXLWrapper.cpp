/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

// The following copyright applies to DecodeJpegXlOneShot() and
// EncodeJpegXLOneShot(): (note that the original source code has been
// modified for Siri's purposes).

// Copyright (c) the JPEG XL Project Authors. All rights reserved.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE_JPEGXL file.

// This C++ example decodes a JPEG XL image in one shot (all input bytes
// available at once). The example outputs the pixels and color information to a
// floating point image and an ICC profile on disk.
//#ifdef HAVE_LIBJXL

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_LIBJXL

#include <inttypes.h>
#include <jxl/decode.h>
#include <jxl/decode_cxx.h>
#include <jxl/encode.h>
#include <jxl/encode_cxx.h>
#include <jxl/thread_parallel_runner.h>
#include <jxl/thread_parallel_runner_cxx.h>
#include <jxl/resizable_parallel_runner.h>
#include <jxl/resizable_parallel_runner_cxx.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <vector>

/** Decodes JPEG XL image to floating point pixels and ICC Profile. Pixel are
 * stored as floating point, as interleaved RGB (3 floating point values per
 * pixel), line per line from top to bottom.  Pixel values have nominal range
 * 0..1 but may go beyond this range for HDR or wide gamut. The ICC profile
 * describes the color format of the pixel data.
 */
bool DecodeJpegXlOneShot(const uint8_t* jxl, size_t size,
                         std::vector<float>* pixels, size_t* xsize,
                         size_t* ysize, size_t* zsize, uint8_t* bitdepth, std::vector<uint8_t>* icc_profile) {
  // Multi-threaded parallel runner.
  auto runner = JxlResizableParallelRunnerMake(nullptr);

  auto dec = JxlDecoderMake(nullptr);
  if (JXL_DEC_SUCCESS !=
      JxlDecoderSubscribeEvents(dec.get(), JXL_DEC_BASIC_INFO |
                                               JXL_DEC_COLOR_ENCODING |
                                               JXL_DEC_FULL_IMAGE)) {
    fprintf(stderr, "JxlDecoderSubscribeEvents failed\n");
    return false;
  }

  if (JXL_DEC_SUCCESS != JxlDecoderSetParallelRunner(dec.get(),
                                                     JxlResizableParallelRunner,
                                                     runner.get())) {
    fprintf(stderr, "JxlDecoderSetParallelRunner failed\n");
    return false;
  }

  JxlBasicInfo info;

  JxlDecoderSetInput(dec.get(), jxl, size);
  JxlDecoderCloseInput(dec.get());
  JxlPixelFormat format = {3, JXL_TYPE_FLOAT, JXL_NATIVE_ENDIAN, 0};

  for (;;) {
    JxlDecoderStatus status = JxlDecoderProcessInput(dec.get());

    if (status == JXL_DEC_ERROR) {
      fprintf(stderr, "Decoder error\n");
      return false;
    } else if (status == JXL_DEC_NEED_MORE_INPUT) {
      fprintf(stderr, "Error, already provided all input\n");
      return false;
    } else if (status == JXL_DEC_BASIC_INFO) {
      if (JXL_DEC_SUCCESS != JxlDecoderGetBasicInfo(dec.get(), &info)) {
        fprintf(stderr, "JxlDecoderGetBasicInfo failed\n");
        return false;
      }
      format = {info.num_color_channels, JXL_TYPE_FLOAT, JXL_NATIVE_ENDIAN, 0};
      *xsize = info.xsize;
      *ysize = info.ysize;
      *zsize = info.num_color_channels;
      *bitdepth = info.bits_per_sample;
      JxlResizableParallelRunnerSetThreads(
          runner.get(),
          JxlResizableParallelRunnerSuggestThreads(info.xsize, info.ysize));
    } else if (status == JXL_DEC_COLOR_ENCODING) {
      // Get the ICC color profile of the pixel data
      size_t icc_size;
      if (JXL_DEC_SUCCESS !=
          JxlDecoderGetICCProfileSize(dec.get(), &format, JXL_COLOR_PROFILE_TARGET_DATA,
                                      &icc_size)) {
        fprintf(stderr, "JxlDecoderGetICCProfileSize failed\n");
        return false;
      }
      icc_profile->resize(icc_size);
      if (JXL_DEC_SUCCESS != JxlDecoderGetColorAsICCProfile(
                                 dec.get(), &format, JXL_COLOR_PROFILE_TARGET_DATA,
                                 icc_profile->data(), icc_profile->size())) {
        fprintf(stderr, "JxlDecoderGetColorAsICCProfile failed\n");
        return false;
      }
    } else if (status == JXL_DEC_NEED_IMAGE_OUT_BUFFER) {
      size_t buffer_size;
      if (JXL_DEC_SUCCESS !=
          JxlDecoderImageOutBufferSize(dec.get(), &format, &buffer_size)) {
        fprintf(stderr, "JxlDecoderImageOutBufferSize failed\n");
        return false;
      }
      if (buffer_size != *xsize * *ysize * 12) {
        fprintf(stderr, "Invalid out buffer size %" PRIu64 " %" PRIu64 "\n",
                static_cast<uint64_t>(buffer_size),
                static_cast<uint64_t>(*xsize * *ysize * 12));
        return false;
      }
      pixels->resize(*xsize * *ysize * *zsize);
      void* pixels_buffer = (void*)pixels->data();
      size_t pixels_buffer_size = pixels->size() * sizeof(float);
      if (JXL_DEC_SUCCESS != JxlDecoderSetImageOutBuffer(dec.get(), &format,
                                                         pixels_buffer,
                                                         pixels_buffer_size)) {
        fprintf(stderr, "JxlDecoderSetImageOutBuffer failed\n");
        return false;
      }
    } else if (status == JXL_DEC_FULL_IMAGE) {
      // Nothing to do. Do not yet return. If the image is an animation, more
      // full frames may be decoded. This example only keeps the last one.
    } else if (status == JXL_DEC_SUCCESS) {
      // All decoding successfully finished.
      // It's not required to call JxlDecoderReleaseInput(dec.get()) here since
      // the decoder will be destroyed.
      return true;
    } else {
      fprintf(stderr, "Unknown decoder status\n");
      return false;
    }
  }
}

extern "C" int DecodeJpegXlOneShotWrapper(const uint8_t* jxl, size_t size,
                         float** pixels, size_t* xsize,
                         size_t* ysize, size_t* zsize, uint8_t* bitdepth,
                         uint8_t** icc_profile, size_t *icc_profile_length) {
    std::vector<float> vec_pixels;
    std::vector<uint8_t> vec_icc_profile;
    if (!DecodeJpegXlOneShot(jxl, size,
                         &vec_pixels, xsize,
                         ysize, zsize, bitdepth,
                         &vec_icc_profile))
        return -1;
    float *array = (float*) malloc(vec_pixels.size() * sizeof(float));
    *pixels = array;
    uint8_t *icc_array = (uint8_t*) malloc(vec_icc_profile.size() * sizeof(uint8_t));
    *icc_profile = icc_array;
    memcpy(*pixels, vec_pixels.data(), vec_pixels.size() * sizeof(float));
    memcpy(*icc_profile, vec_icc_profile.data(), vec_icc_profile.size() * sizeof(uint8_t));
    return 0;
}

/**
 * Compresses the provided pixels.
 *
 * @param pixels input pixels
 * @param xsize width of the input image
 * @param ysize height of the input image
 * @param compressed will be populated with the compressed bytes
 */
bool EncodeJxlOneshot(const std::vector<uint8_t>& pixels, const uint32_t xsize,
                      const uint32_t ysize, const uint32_t zsize, const uint8_t bitdepth,
                      std::vector<uint8_t>* compressed, const uint32_t effort, const float distance) {
  auto enc = JxlEncoderMake(/*memory_manager=*/nullptr);
#ifdef HAVE_LIBJXL_THREADS
  auto runner = JxlThreadParallelRunnerMake(
      /*memory_manager=*/nullptr,
      JxlThreadParallelRunnerDefaultNumWorkerThreads());
  if (JXL_ENC_SUCCESS != JxlEncoderSetParallelRunner(enc.get(),
                                                     JxlThreadParallelRunner,
                                                     runner.get())) {
    fprintf(stderr, "JxlEncoderSetParallelRunner failed\n");
    return false;
  }
#endif

  JxlPixelFormat pixel_format;
  switch (bitdepth) {
    case 8:
      pixel_format = {zsize, JXL_TYPE_UINT8, JXL_NATIVE_ENDIAN, 0};
      break;
    case 16:
      pixel_format = {zsize, JXL_TYPE_UINT16, JXL_NATIVE_ENDIAN, 0};
      break;
    case 32:
      pixel_format = {zsize, JXL_TYPE_FLOAT, JXL_NATIVE_ENDIAN, 0};
  }

  JxlBasicInfo basic_info;
  JxlEncoderInitBasicInfo(&basic_info);
  basic_info.xsize = xsize;
  basic_info.ysize = ysize;
  basic_info.bits_per_sample = bitdepth;
  basic_info.exponent_bits_per_sample = 8;
  basic_info.uses_original_profile = distance == 0.0 ? JXL_TRUE : JXL_FALSE;
  if (JXL_ENC_SUCCESS != JxlEncoderSetBasicInfo(enc.get(), &basic_info)) {
    fprintf(stderr, "JxlEncoderSetBasicInfo failed\n");
    return false;
  }

  JxlColorEncoding color_encoding = {};
  JxlColorEncodingSetToSRGB(&color_encoding,
                            /*is_gray=*/pixel_format.num_channels < 3);
  if (JXL_ENC_SUCCESS !=
      JxlEncoderSetColorEncoding(enc.get(), &color_encoding)) {
    fprintf(stderr, "JxlEncoderSetColorEncoding failed\n");
    return false;
  }

  JxlEncoderFrameSettings* frame_settings =
      JxlEncoderFrameSettingsCreate(enc.get(), nullptr);

  if (effort > 0 && effort < 10)
    JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_EFFORT,
                                  effort);
  else
      fprintf(stderr, "JxlEncoderFrameSettings: effort is outside limits. Using default (= 7)\n");

  if (JXL_ENC_SUCCESS != JxlEncoderSetFrameDistance(frame_settings, distance))
    fprintf(stderr, "Error setting permissible distance: continuing with default.\n");

  if (JXL_ENC_SUCCESS !=
      JxlEncoderAddImageFrame(frame_settings, &pixel_format,
                              (void*)pixels.data(),
                              pixels.size())) {
    fprintf(stderr, "JxlEncoderAddImageFrame failed\n");
    return false;
  }
  JxlEncoderCloseInput(enc.get());

  compressed->resize(64);
  uint8_t* next_out = compressed->data();
  size_t avail_out = compressed->size() - (next_out - compressed->data());
  JxlEncoderStatus process_result = JXL_ENC_NEED_MORE_OUTPUT;
  while (process_result == JXL_ENC_NEED_MORE_OUTPUT) {
    process_result = JxlEncoderProcessOutput(enc.get(), &next_out, &avail_out);
    if (process_result == JXL_ENC_NEED_MORE_OUTPUT) {
      size_t offset = next_out - compressed->data();
      compressed->resize(compressed->size() * 2);
      next_out = compressed->data() + offset;
      avail_out = compressed->size() - offset;
    }
  }
  compressed->resize(next_out - compressed->data());
  if (JXL_ENC_SUCCESS != process_result) {
    fprintf(stderr, "JxlEncoderProcessOutput failed\n");
    return false;
  }

  return true;
}

extern "C" int EncodeJpegXlOneshotWrapper(const uint8_t* pixels, const uint32_t xsize,
                      const uint32_t ysize, const uint32_t zsize, const uint8_t bitdepth,
                      void** compressed, size_t* compressed_length, uint32_t effort,
                      const double distance) {
    int datasize = bitdepth / 8;
    std::vector<uint8_t> vec_pixels(pixels, pixels + xsize * ysize * zsize * datasize);
    std::vector<uint8_t> vec_compressed;

    int retval = (!EncodeJxlOneshot(vec_pixels, xsize,
                      ysize, zsize, bitdepth,
                      &vec_compressed, effort, distance)) ? 1 : 0;
    void* array = (void*) malloc(vec_compressed.size() * datasize);
    *compressed = array;
    memcpy(*compressed, vec_compressed.data(), vec_compressed.size() * sizeof(uint8_t));
    *compressed_length = vec_compressed.size();
     return retval;
}

#endif
