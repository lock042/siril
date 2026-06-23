/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
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

/*
 * Fast first-frame extractor for video files (AVI / MP4 / MOV / ...).
 *
 * The Siril sequence pipeline opens video files through FFMS2, which
 * builds a full-file index up-front to provide random-access seeking
 * during stacking.  Indexing a multi-GB AVI takes several seconds —
 * fine when you're about to stack thousands of frames, terrible when
 * you just clicked the file in the browser to see what it looks like.
 *
 * This extractor uses libavformat / libavcodec directly to decode just
 * the first frame: open the container, find the video stream, read
 * packets until one decodes into a frame, convert to RGB24, return.
 * No indexing, no full-file scan — typically just the first few hundred
 * kilobytes are read.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <glib.h>

#ifdef HAVE_FFMPEG

#include <libavformat/avformat.h>
#include <libavcodec/avcodec.h>
#include <libswscale/swscale.h>
#include <libavutil/pixdesc.h>

#include "io/avi_preview.h"

guchar *extract_thumbnail_from_avi(const char *filename, gchar **descr,
                                    int *width_out, int *height_out) {
	if (descr)       *descr = NULL;
	if (width_out)   *width_out = 0;
	if (height_out)  *height_out = 0;

	AVFormatContext *fmt = NULL;
	if (avformat_open_input(&fmt, filename, NULL, NULL) < 0)
		return NULL;
	if (avformat_find_stream_info(fmt, NULL) < 0) {
		avformat_close_input(&fmt);
		return NULL;
	}

	/* Locate the first video stream. */
	int vstream_idx = -1;
	AVCodecParameters *codecpar = NULL;
	for (unsigned int i = 0; i < fmt->nb_streams; i++) {
		if (fmt->streams[i]->codecpar->codec_type == AVMEDIA_TYPE_VIDEO) {
			vstream_idx = (int) i;
			codecpar = fmt->streams[i]->codecpar;
			break;
		}
	}
	if (vstream_idx < 0 || !codecpar) {
		avformat_close_input(&fmt);
		return NULL;
	}

	const AVCodec *codec = avcodec_find_decoder(codecpar->codec_id);
	if (!codec) {
		avformat_close_input(&fmt);
		return NULL;
	}
	AVCodecContext *cctx = avcodec_alloc_context3(codec);
	if (!cctx) {
		avformat_close_input(&fmt);
		return NULL;
	}
	if (avcodec_parameters_to_context(cctx, codecpar) < 0
	    || avcodec_open2(cctx, codec, NULL) < 0) {
		avcodec_free_context(&cctx);
		avformat_close_input(&fmt);
		return NULL;
	}

	AVPacket *pkt = av_packet_alloc();
	AVFrame  *frame = av_frame_alloc();
	if (!pkt || !frame) {
		if (pkt) av_packet_free(&pkt);
		if (frame) av_frame_free(&frame);
		avcodec_free_context(&cctx);
		avformat_close_input(&fmt);
		return NULL;
	}

	/* Pull packets off the wire until the decoder hands us a frame.
	 * A handful of early packets may be I-frame headers or B-frame
	 * deltas that need more context — that's why we loop rather than
	 * just decoding the first packet we get. */
	gboolean have_frame = FALSE;
	while (av_read_frame(fmt, pkt) >= 0) {
		if (pkt->stream_index != vstream_idx) {
			av_packet_unref(pkt);
			continue;
		}
		int send_r = avcodec_send_packet(cctx, pkt);
		av_packet_unref(pkt);
		if (send_r < 0) break;
		int recv_r = avcodec_receive_frame(cctx, frame);
		if (recv_r == 0) { have_frame = TRUE; break; }
		if (recv_r == AVERROR(EAGAIN)) continue;
		if (recv_r < 0) break;
	}
	/* Flush any decoder-buffered frames if we exited the read loop
	 * without getting one (the file may be very short). */
	if (!have_frame) {
		avcodec_send_packet(cctx, NULL);
		if (avcodec_receive_frame(cctx, frame) == 0)
			have_frame = TRUE;
	}

	guchar *rgb = NULL;
	if (have_frame && frame->width > 0 && frame->height > 0) {
		int w = frame->width;
		int h = frame->height;
		struct SwsContext *sws = sws_getContext(w, h, (enum AVPixelFormat) frame->format,
			w, h, AV_PIX_FMT_RGB24,
			SWS_BILINEAR, NULL, NULL, NULL);
		if (sws) {
			rgb = malloc((size_t) w * h * 3);
			if (rgb) {
				uint8_t *dst[1]   = { rgb };
				int      stride[1] = { w * 3 };
				sws_scale(sws, (const uint8_t * const *) frame->data,
				          frame->linesize, 0, h, dst, stride);
				if (width_out)  *width_out  = w;
				if (height_out) *height_out = h;
				if (descr) {
					/* Translation-friendly description.  We use ASCII
					 * substitutes — the file_browser.c default-preview
					 * already handles type / size separately, so this
					 * string just needs to carry dimensions and frame
					 * count back to the caller. */
					int nb_frames = (int) fmt->streams[vstream_idx]->nb_frames;
					if (nb_frames > 0) {
						*descr = g_strdup_printf("%d x %d pixels\n%d frames",
							w, h, nb_frames);
					} else {
						*descr = g_strdup_printf("%d x %d pixels", w, h);
					}
				}
			}
			sws_freeContext(sws);
		}
	}

	av_frame_free(&frame);
	av_packet_free(&pkt);
	avcodec_free_context(&cctx);
	avformat_close_input(&fmt);
	return rgb;
}

#else  /* !HAVE_FFMPEG */

#include "io/avi_preview.h"

guchar *extract_thumbnail_from_avi(const char *filename, gchar **descr,
                                    int *width_out, int *height_out) {
	(void) filename;
	if (descr)      *descr = NULL;
	if (width_out)  *width_out = 0;
	if (height_out) *height_out = 0;
	return NULL;
}

#endif /* HAVE_FFMPEG */
