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
    """
    Represents a single data series for plotting.

    Members:
        x_coords: Either a List[float] or a np.ndarray containing the values
        for the x coordinates for this series

        y_coords: Either a List[float] or a np.ndarray containing the values
        for the y coordinates for this series

        label: A str containing a label for the series (shown in the plot
        legend)

        plot_type: a PlotType setting the type of marks to use

        n_error: Either a List[float] or a np.ndarray containing values for
        the y-axis negative errors for this series

        p_error: Either a List[float] or a np.ndarray containing values for
        the y-axis positive errors for this series
    """

    def __init__(
        self,
        x_coords: Union[List[float], np.ndarray],
        y_coords: Union[List[float], np.ndarray],
        label: Optional[str] = None,
        plot_type: Optional[PlotType] = PlotType.LINES,
        n_error: Optional[Union[List[float], np.ndarray]] = None,
        p_error: Optional[Union[List[float], np.ndarray]] = None
    ):
        """
        Represents a single data series for plotting.

        Args:
            x_coords: X-coordinates of the data series
            y_coords: Y-coordinates of the data series
            label: Label for the series (optional)
            plot_type: Type of plot for this series (optional, default is LINES)
            n_error: Y-axis negative error for error bars (optional)
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
        self.n_error = np.asarray(n_error, dtype=np.float32) if n_error is not None else None
        self.p_error = np.asarray(p_error, dtype=np.float32) if p_error is not None else None

        # Validate error bar lengths if provided
        if self.n_error is not None and len(self.n_error) != len(self.y_coords):
            raise ValueError("n_error must have the same length as y_coords")
        if self.p_error is not None and len(self.p_error) != len(self.y_coords):
            raise ValueError("p_error must have the same length as y_coords")

class PlotData:
    """
    Metadata container for plot configuration. The actual series data are
    held in SeriesData objects and can be added using the Class methods
    add_series or add_series_obj after initialization of the PlotData.

    Members:
        title: Plot title

        xlabel: X-axis label

        ylabel: Y-axis label

        savename: Save filename (extension is added automatically)

        show_legend: bool indicating whether to show legend

        datamin: List [xmin, ymin] forcing the bottom left coordinate to show.
        If omitted, the range is set to the data range.

        datamax: List [xmax, ymax] forcing the top right coordinate to show.
        If omitted, the range is set to the data range.
    """

    def __init__(
        self,
        title: Optional[str] = "Data Plot",
        xlabel: Optional[str] = "X",
        ylabel: Optional[str] = "Y",
        savename: Optional[str] = "plot",
        show_legend: Optional[bool] = True,
        datamin: Optional[List[float]] = None,
        datamax: Optional[List[float]] = None
    ):
        """
        Metadata container for plot configuration.

        Args:
            title: Plot title
            xlabel: X-axis label
            ylabel: Y-axis label
            savename: Save filename (extension is added automatically)
            show_legend: bool indicating whether to show legend
            datamin: List [xmin, ymin] forcing the bottom left coordinate to show
            datamax: List [xmax, ymax] forcing the top right coordinate to show
        """
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.savename = savename
        self.show_legend = show_legend
        self.series_data: List[SeriesData] = []
        self.datamin = datamin
        self.datamax = datamax

        # Validate datamin and datamax
        if self.datamin is not None:
            if not isinstance(self.datamin, list):
                raise TypeError("datamin must be a list of 2 numeric values (integers or floats)")
            if len(self.datamin) != 2:
                raise ValueError("datamin must contain exactly 2 numeric values (integers or floats)")
            if not all(isinstance(x, (int, float)) for x in self.datamin):
                raise TypeError("datamin must contain only numeric values (integers or floats)")

        if self.datamax is not None:
            if not isinstance(self.datamax, list):
                raise TypeError("datamax must be a list of 2 numeric values (integers or floats)")
            if len(self.datamax) != 2:
                raise ValueError("datamax must contain exactly 2 numeric values (integers or floats)")
            if not all(isinstance(x, (int, float)) for x in self.datamax):
                raise TypeError("datamax must contain only numeric values (integers or floats)")

    def add_series(
        self,
        x_coords: Union[List[float], np.ndarray],
        y_coords: Union[List[float], np.ndarray],
        label: Optional[str] = None,
        plot_type: Optional[PlotType] = PlotType.LINES,
        n_error: Optional[Union[List[float], np.ndarray]] = None,
        p_error: Optional[Union[List[float], np.ndarray]] = None
    ) -> SeriesData:
        """
        Add a new series to the plot metadata.

        Returns:
            SeriesData: the created SeriesData object for further manipulation if needed.
        """
        series = SeriesData(
            x_coords,
            y_coords,
            label,
            plot_type,
            n_error,
            p_error
        )
        self.series_data.append(series)
        return series

    def add_series_obj(self, series: SeriesData) -> None:
        """
        Add a pre-created SeriesData object to the plot metadata.

        Returns: None
        """
        self.series_data.append(series)

class _PlotSerializer:
    @staticmethod
    def _serialize_plot_data(Plot_data: PlotData) -> Tuple[bytes, int]:
        """
        Serialize plot data for shared memory transfer using network byte order.

        Args:
            Plot_data: PlotData object containing plot configuration

        Returns:
            Tuple of serialized bytes and total length
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

        # If datamin is set, serialize it
        serialized += struct.pack('!?', Plot_data.datamin is not None)
        if Plot_data.datamin is not None:
            serialized += struct.pack('!dd', Plot_data.datamin[0], Plot_data.datamin[1])

        # If datamax is set, serialize it
        serialized += struct.pack('!?', Plot_data.datamax is not None)
        if Plot_data.datamax is not None:
            serialized += struct.pack('!dd', Plot_data.datamax[0], Plot_data.datamax[1])

        for series in Plot_data.series_data:
            with_errors = series.n_error is not None or series.p_error is not None
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
                    if series.n_error is not None and i < len(series.n_error):
                        serialized += struct.pack('!d', series.n_error[i])
                    else:
                        serialized += struct.pack('!d', 0.0)

                    # Serialize positive error (if exists, otherwise 0)
                    if series.p_error is not None and i < len(series.p_error):
                        serialized += struct.pack('!d', series.p_error[i])
                    else:
                        serialized += struct.pack('!d', 0.0)

        return serialized, len(serialized)
