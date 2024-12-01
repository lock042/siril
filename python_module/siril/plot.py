# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np
import struct
import sys
import os
import time
from typing import Union, Optional, List, Tuple
from enum import IntEnum

class PlotType(IntEnum):
    POINTS = 0
    MARKS = 1
    HYPHENS = 2
    LINES = 3
    LINESPOINTS = 4
    LINESMARKS = 5
    LINESHYPHENS= 6

class SeriesData:
    def __init__(
        self,
        x_coords: Union[List[float], np.ndarray],
        y_coords: Union[List[float], np.ndarray],
        label: Optional[str] = None,
        plot_type: Optional[PlotType] = PlotType.LINES,
        m_error: Optional[Union[List[float], np.ndarray]] = None,
        p_error: Optional[Union[List[float], np.ndarray]] = None
    ):
        """
        Represents a single data series for plotting.

        Args:
            x_coords: X-coordinates of the data series
            y_coords: Y-coordinates of the data series
            label: Label for the series (optional)
            plot_type: Type of plot for this series (default: LINES)
            m_error: Y-axis negative error for error bars (optional)
            p_error: Y-axis positive error for error bars (optional)
        """
        # Convert inputs to numpy arrays
        self.x_coords = np.asarray(x_coords, dtype=np.float32)
        self.y_coords = np.asarray(y_coords, dtype=np.float32)

        # Validate coordinate lengths
        if len(self.x_coords) != len(self.y_coords):
            raise ValueError("x and y coordinates must have the same length")

        # Set optional parameters
        self.label = label or "Data Series"
        self.plot_type = plot_type

        # Handle error bars
        self.m_error = np.asarray(m_error, dtype=np.float32) if m_error is not None else None
        self.p_error = np.asarray(p_error, dtype=np.float32) if p_error is not None else None

        # Validate error bar lengths if provided
        if self.m_error is not None and len(self.m_error) != len(self.y_coords):
            raise ValueError("m_error must have the same length as y_coords")
        if self.p_error is not None and len(self.p_error) != len(self.y_coords):
            raise ValueError("p_error must have the same length as y_coords")

class PlotData:
    def __init__(
        self,
        title: Optional[str] = "Data Plot",
        xlabel: Optional[str] = "X",
        ylabel: Optional[str] = "Y",
        savename: Optional[str] = "plot",
        show_legend: Optional[bool] = True
    ):
        """
        Metadata container for plot configuration.

        Args:
            title: Plot title
            xlabel: X-axis label
            ylabel: Y-axis label
            savename: Save filename
            show_legend: Boolean indicating whether to show legend
        """
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.savename = savename
        self.show_legend = show_legend
        self.series_data: List[SeriesData] = []

    def add_series(
        self,
        x_coords: Union[List[float], np.ndarray],
        y_coords: Union[List[float], np.ndarray],
        label: Optional[str] = None,
        plot_type: Optional[PlotType] = PlotType.LINES,
        m_error: Optional[Union[List[float], np.ndarray]] = None,
        p_error: Optional[Union[List[float], np.ndarray]] = None
    ) -> SeriesData:
        """
        Add a new series to the plot metadata.

        Returns the created SeriesData object for further manipulation if needed.
        """
        series = SeriesData(
            x_coords,
            y_coords,
            label,
            plot_type,
            m_error,
            p_error
        )
        self.series_data.append(series)
        return series

    def add_series_obj(self, series: SeriesData) -> None:
        """
        Add a pre-created SeriesData object to the plot metadata.
        """
        self.series_data.append(series)

class PlotSerializer:
    @staticmethod
    def _serialize_plot_data(Plot_data: PlotData) -> Tuple[bytes, int]:
        """
        Serialize plot data for shared memory transfer using network byte order.
        Args:
            Plot_data: PlotData object containing plot configuration
        Returns: Tuple of serialized bytes and total length
        """
        def encode_null_string(s):
            return s.encode('utf-8') + b'\x00'

        serialized = b''
        serialized += encode_null_string(Plot_data.title or "")
        serialized += encode_null_string(Plot_data.xlabel or "")
        serialized += encode_null_string(Plot_data.ylabel or "")
        serialized += encode_null_string(Plot_data.savename or "")

        # Pack boolean and number of series
        serialized += struct.pack('!?', Plot_data.show_legend)
        serialized += struct.pack('!I', len(Plot_data.series_data))

        for series in Plot_data.series_data:
            with_errors = series.m_error is not None or series.p_error is not None
            serialized += encode_null_string(series.label)
            serialized += struct.pack('!?', with_errors)
            serialized += struct.pack('!I', len(series.x_coords))
            serialized += struct.pack('!I', series.plot_type.value)

            # Serialize coordinates with optional error bars for each point
            for i in range(len(series.x_coords)):
                # Serialize x and y coordinates
                serialized += struct.pack('!dd', series.x_coords[i], series.y_coords[i])

                if with_errors:
                    # Serialize negative error (if exists, otherwise 0)
                    if series.m_error is not None and i < len(series.m_error):
                        serialized += struct.pack('!d', series.m_error[i])
                    else:
                        serialized += struct.pack('!d', 0.0)

                    # Serialize positive error (if exists, otherwise 0)
                    if series.p_error is not None and i < len(series.p_error):
                        serialized += struct.pack('!d', series.p_error[i])
                    else:
                        serialized += struct.pack('!d', 0.0)

        return serialized, len(serialized)
