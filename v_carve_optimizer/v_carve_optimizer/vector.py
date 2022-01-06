#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# pylint: disable=unexpected-keyword-arg,E1101,E1130,E0602,R0913,R0201,C0103
#
#
#       Copyright 2009-2014 Pelayo González <pelayogmp@gmail.com>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.


r"""
:mod:`geo.vector` -- Vector and Point arithmetic
================================================

..  module:: geo.vector
    :synopsis: Vector and Point arithmetic and utilities.
..  moduleauthor:: Pelayo González <pelayogmp@gmail.com>
..  currentmodule:: primfitter.geo.vector

..  Numpy documentation urls

..  _ndarray: http://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html
..  _array: http://docs.scipy.org/doc/numpy/reference/generated/numpy.array.html#numpy-array
..  _broadcasting: http://docs.scipy.org/doc/numpy/user/basics.broadcasting.html
..  _matrix: http://docs.scipy.org/doc/numpy/reference/arrays.classes.html#matrix-objects
..  _view: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.view.html#numpy-ndarray-view
..  _dtype: http://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html#data-type-objects-dtype
..  _orthogonal: http://en.wikipedia.org/wiki/Orthonormal_matrix
..  _mlab: http://code.enthought.com/projects/mayavi/documentation.php
..  _Mayavi2: http://code.enthought.com/projects/mayavi/

Utility module for vectors and points arthmetic.

External dependencies
---------------------

..  digraph:: "geo_vector_depend"

        rankdir="LR"

        node [fontname="Helvetica", fontsize="10", shape="box",
              height="0.25", width="0.25", style="setlinewidth(0.5)"]

        edge [fontname="Helvetica", fontsize="10", arrowsize="0.5"
              style="setlinewidth(0.5)", weight="2.0"]

        /* Dependency packages */
        {
            node [fontcolor="#000060"]
            sys [href="http://docs.python.org/library/sys.html"]
            numpy [href="http://docs.scipy.org/doc/numpy/reference/"]
        }

        sys -> vector
        numpy -> vector

Future development
------------------

    * Implement general transformations.
      See http://en.wikipedia.org/wiki/Quaternion

Public functions and classes
----------------------------


"""

from __future__ import (print_function,
                        division,
                        absolute_import,
                        unicode_literals)
import contextlib
import sys
import numpy


@contextlib.contextmanager
def printoptions(*args, **kwargs):
    r"""
    Generator to configure numpyprintoptions and return back to original
    """
    original = numpy.get_printoptions()
    numpy.set_printoptions(*args, **kwargs)
    yield
    numpy.set_printoptions(**original)

def one_dim_vector(vec):
    """
    Ensures that vec is ndim = 1 and shape = (3,)
    """
    _vec = Vector(vec)
    _shape = _vec.shape
    _condition = (_shape == (3,)) or (_shape == (1, 3))
    assert _condition, "In one_dim_vector invalid shape %s" % str(_shape)
    return _vec.reshape(3,)

def check_orthogonal(axis):
    r"""
    Checks if a matrix is orthogonal.

    An orthogonal_ matrix is:

        * Two dimensional.
        * Real.
        * Square.
        * Invertible.
        * Its inverse and its transpose are equal.

    **Arguments**

    ``axis`` *array shape(3, 3)*

        when normalized Must be a orthogonal_ matrix.

    **Returns**

    **matrix shape(3, 3)**

        Orthogonal matrix.

    **Raises**

    ``ValueError``

        If the matrix is not orthogonal.

    ..  seealso::

        numpy matrix_ class.

    """
    axis = Vector(axis).unit
    _is_orthogonal = (axis.ndim == 2) and \
                     (axis.shape[0] == axis.shape[1]) and \
                     numpy.isreal(axis).all()

    if _is_orthogonal:
        axis_mat = numpy.asmatrix(axis)
        # a.T == a.I only for orthogonal matrix
        _is_orthogonal = numpy.allclose(axis_mat.T, axis_mat.I)

    if _is_orthogonal:
        return axis

    raise ValueError("Matrix is not orthogonal.")

class Vector(numpy.ndarray):
    r"""
    Class for 3D vectors and points.

    ..  inheritance-diagram:: Vector

    :class:`Vector` is subclass of numpy ndarray_ [#]_ [#]_, a *vector* is an
    ndarray_ with at most two dimensions, every row in the array defines a
    single position or direction. In :math:`\mathbb{R}^{3}` arrays
    are of *shape(3)* for single vectors or *shape(N, 3)* for *N* vectors.

    All methods are designed to support numpy broadcasting_. So
    operations accepts the shape combinations allowed by numpy broadcasting_
    rules. If rules are not obeyed a ``ValuError`` exception is raised::

        shape mismatch: objects cannot be broadcast to a single shape.

    ..  rubric:: Footnotes

    ..  [#] http://docs.scipy.org/doc/numpy/reference/arrays.classes.html
    ..  [#] http://docs.scipy.org/doc/numpy/user/basics.subclassing.html

    """

    __array_priority__ = 10.0

    def __new__(cls, data, dtype=None, copy=False):
        r"""
        **Example**

        >>> beast = numpy.random.rand(1000000,3)
        >>> %timeit vector.Vector(beast)
        100000 loops, best of 3: 3.52 us per loop
        """
        # Make sure we are working with an array, and copy the data
        # if requested
        obj = numpy.array(data, dtype=dtype, copy=copy)

        # Transform 'subarr' from an ndarray to our new subclass.
        obj = obj.view(cls)

        # Finally, we must return the newly created object:
        return obj

    def __init__(self, data, dtype=None, copy=False):
        r"""

        **Arguments**

        ``data``
            An array definition suitable for numpy array_ function.

        **Keyword arguments**

        ``dtype`` *numpy.dtype* ``(default None)``
            Any valid numpy dtype_.

        ``copy`` *bool* ``(default False)``
            Whether return a copy or a view_ of *data*.

        """

        pass

    def __array_finalize__(self, obj):
        r"""
        Check if compatible.
        """

        if obj is None:
            return

        if obj.ndim > 2:
            raise ValueError("Invalid dimension %d" % obj.ndim)


    def dot(self, vec):
        r"""
        Returns the dot product.
        """
        return (self * Vector(vec)).sum(axis=-1)

    def cross(self, vec):
        r"""
        Return the cross product.
        """
        return Vector(numpy.cross(self, Vector(vec)))

    def triple(self, vec_a, vec_b):
        r"""
        Returns the triple product.

        The triple product is defined as::

            (self.cross(vec_a)).dot(vec_b)

        """
        return (self.cross(Vector(vec_a))).dot(Vector(vec_b))

    def angle(self, vec):
        r"""
        Returns the angle between *self* and *vec* in radians.
        """

        vec = Vector(vec)
        mod_prod = self.leng * vec.leng
        return numpy.arccos(self.dot(vec)/mod_prod)

    def degrees(self, vec):
        r"""
        Returns the angle between *self* and *vec* in degrees.
        """
        return numpy.rad2deg(self.angle(vec))

    def orthogonal(self, vec):
        r"""
        Returns the projection of *vec* onto the plane normal to *self*
        containing the origin.

        So if :math:`(\vec{V})` is *vec* and :math:`(\vec{V_{u}})` is *self.unit*:

        .. math::

            \vec{V_{p}} = \vec{V} - (\vec{V} \bullet \vec{V_{u}}) \vec{V_{u}}

        Where :math:`\bullet` denotes dot product.
        """
        vec = Vector(vec)
        self_u = self.unit
        dot_p = self_u.dot(vec)
        if self.ndim != 1 or vec.ndim != 1:
            dot_p = numpy.tile(dot_p, (3, 1)).T

        return vec - dot_p * self_u

    def project(self, other):
        r"""
        Returns the projection of *self* onto *other* entity.

        ..  note::

            *other* class must implement the method ``project_point()``.

        """

        try:
            other = Vector(other, copy=False)
        except TypeError:
            pass
        if isinstance(other, Vector):
            other_u = other.unit
            dot_p = self.dot(other_u)
            if self.ndim != 1 or other_u.ndim != 1:
                dot_p = numpy.tile(dot_p, (3, 1)).T
            return dot_p * other_u
        elif hasattr(other, project_point):
            return other.project_point(self)
        else:
            raise NotImplementedError("%s does not implement project_point method" \
                % other.__class__)

    def distance(self, other):
        r"""
        Returns the distance to *other* entity.

        ..  note::

            *other* class must implement the method ``distance_to_point()``.

        """
        try:
            other = Vector(other, copy=True)
        except TypeError:
            pass
        if isinstance(other, Vector):
            return (self - other).leng
        elif hasattr(other, distance_to_point):
            return other.distance_to_point(self)
        else:
            raise NotImplementedError( \
                "%s does not implement distance_to_point method" \
                % other.__class__)

    def transform(self, axes, origin=None, reverse=False, check_axes=True, copy=True):
        r"""
        Returns a copy transformed by *axes* and *origin*.

        The transformed :class:`Vector` coordinates are measured from
        the reference defined by *axes* and *origin*. But the opposite
        behavior is obtained with the *reverse* flag.

        **Arguments**

        ``axes`` *array shape(3, 3)*
            Vectors of transformation by rows::

                i = axes[0], j = axes[1], k = axes[2]

            Should be a valid orthogonal_ *(orthonormal)* matrix, but see
            ``check_axes`` and :func:`check_orthogonal`.

        **Keywords arguments**

        ``origin`` *array shape(3, 3)* ``(default None)``
            Origin of the transformation, if ``None`` only rotation
            is applied.

        ``reverse`` *bool* ``(default False)``
            Unapplies the transform.

        ``check_axes`` *bool* ``(default True)``
            Check if ``axes`` is a valid orthogonal_ matrix.

        ``copy`` *bool* ``(default True)``
            Copy the transformed vector in another vector.

        **Returns**

        *Vector same shape as self*
            Transformed points. May return ``None`` if ``axes`` is not
            orthogonal and ``check_axes`` is ``True``.

        ..  seealso::

            :func:`check_orthogonal`.

        """
        if check_axes:
            try:
                axes = check_orthogonal(axes)
            except ValueError:
                print("Vector.transform: Invalid transformation.",
                      file=sys.stderr)
                return None

        axes = numpy.asmatrix(axes)
        if origin is not None:
            origin = Vector(origin)
        if reverse:
            if not origin is None:
                origin = -origin.transform(axes, check_axes=False)
        else:
            axes = axes.T
        # Transformation are easy but tricky.
        if origin is None:
            transformed = Vector(numpy.matrix(self) * axes)
        else:
            transformed = Vector(numpy.matrix(self - origin) * axes)

        if not copy:
            self = transformed

        return transformed

    def scale(self, factor, yfactor=None, zfactor=None, origin=None, copy=False):
        r"""
        Scale a vector.
        """
        if (yfactor is None) and (zfactor is None):
            if origin is None:
                origin = [0.0, 0.0, 0.0]

            origin = Vector(origin)
            scaled = (self - origin) * factor + origin
            if not copy:
                self = scaled
            return scaled

        else:
            raise NotImplementedError("Not implemented in this version")

    def set_az_el(self, az, el, deg=False):
        """
        Sets the vector az el angles

        az, el may be single values or arrays of shape self.shape[0]
        """
        if deg:
            az = numpy.deg2rad(az)
            el = numpy.deg2rad(el)

        self.z = numpy.sin(el)
        _proj = numpy.cos(el)
        self.x = _proj * numpy.cos(az)
        self.y = _proj * numpy.sin(az)

    def _get_cyl(self, deg=False):
        r"""
        Access to cylindrical coordinates as numpy.array.

        Cylindrical coor. = ``[az, rad, z]``, azimuth in radians.
        """

        if self.ndim == 1:
            return Vector([self._get_az(deg), self.rad, self.z])
        return Vector(numpy.vstack(\
                          [self._get_az(deg), self.rad, self.z]).T)

    def _get_cyl_deg(self):
        r"""
        Access to cylindrical coordinates as numpy.array.

        Cylindrical coor. = ``[az, rad, z]``, azimuth in degrees.
        """
        return self._get_cyl(deg=True)

    def _get_sp(self, deg=False):
        r"""
        Access to spherical coordinates as numpy.array.

        Spherical coor. = ``[az, z_angle, leng]``, azimuth and elevation in
        radians.

        """

        if self.ndim == 1:
            return Vector([self._get_az(deg), self._get_z_angle(deg), \
                self.leng])
        return Vector(numpy.vstack( \
                [self._get_az(deg), self._get_z_angle(deg), self.leng]).T)

    def _get_sp_deg(self):
        r"""
        Access to spherical coordinates as numpy.array.

        Spherical coor. = ``[az, z_angle, leng]``, azimuth and elevation in
        degrees.
        """
        return self._get_sp(deg=True)

    def _get_x(self):
        r"""
        Returns the ``X`` coordinates.
        """
        if self.ndim == 1:
            return self[0]
        return self[:, 0]

    def _set_x(self, value):
        if self.ndim == 1:
            self[0] = value
        else:
            self[:, 0] = value

    def _get_y(self):
        r"""
        Returns the ``Y`` coordinates.
        """
        if self.ndim == 1:
            return self[1]
        return self[:, 1]

    def _set_y(self, value):
        if self.ndim == 1:
            self[1] = value
        else:
            self[:, 1] = value

    def _get_z(self):
        r"""
        Returns the ``Z`` coordinates.
        """
        if self.ndim == 1:
            return self[2]
        return self[:, 2]

    def _set_z(self, value):
        if self.ndim == 1:
            self[2] = value
        else:
            self[:, 2] = value

    def _get_leng_sq(self):
        r"""
        Returns the squared lengths
        """
        return self.dot(self)

    def _get_leng(self):
        r"""
        Returns the lengths (spherical radius).
        """
        return numpy.sqrt(self.leng_sq)

    def _get_az(self, deg=False):
        r"""
        Returns the azimuth angle in radians.

        Measured from ``X+`` to the projection onto ``XY`` plane.
        """

        if deg:
            per_rad = 180.0/numpy.pi
        else:
            per_rad = 1.0
        return per_rad * numpy.arctan2(self.y, self.x)


    def _get_az_deg(self):
        r"""
        Returns the azimuth angle in degrees.

        Measured from ``X+`` to the projection onto ``XY`` plane.
        """
        return self._get_az(deg=True)

    def _get_rad(self):
        r"""
        Returns the cylindrical radius.

        The lenght of the projection onto ``XY`` plane.
        """
        if self.ndim == 1:
            return self[:2].leng
        return self[:, :2].leng

    def _get_z_angle(self, deg=False):
        r"""
        Returns the angle with ``Z+`` in radians.
        """
        if deg:
            per_rad = 180.0/numpy.pi
        else:
            per_rad = 1.0

        return per_rad * self.angle([0.0, 0.0, 1.0])

    def _get_z_angle_deg(self):
        r"""
        Returns the angle with ``Z+`` in degrees.

        """
        return self._get_z_angle(deg=True)

    def _get_el(self, deg=False):
        r"""
        Returns the spherical elevation in radians.

        Measured from the projection onto ``XY`` plane to *self*,
        this is the complementary of :meth:`z_angle`.
        """
        if deg:
            per_rad = 180.0/numpy.pi
        else:
            per_rad = 1.0
        _denom = numpy.sqrt(numpy.square(self.x) + numpy.square(self.y))
        return per_rad * numpy.arctan2(self.z, _denom)

    def _get_el_deg(self):
        r"""
        Returns the spherical elevation in degrees.

        """
        return self._get_el(deg=True)

    def _get_unit(self):
        r"""
        Returns a :class:`Vector` of unary vectors.
        """
        if self.ndim == 1:
            return self / self.leng

        return self / numpy.tile(self.leng, (3, 1)).T

    def _get_centroid(self):
        r"""
        Returns the *centroid (center of gravitiy)* of *self*.
        """
        if self.ndim == 1:
            return self

        return self.mean(axis=0)

    def _get_bbox(self):
        r"""
        Returns the bounding box
        """
        _min = [self.x.min(), self.y.min(), self.z.min()]
        _max = [self.x.max(), self.y.max(), self.z.max()]
        return Vector([_min, _max])

    def _get_bbox_center(self):
        r"""
        Returns the center of the bounding box
        """
        _bb = self.bbox
        _center = (_bb[0] + _bb[1])/2.0
        return _center


    def get_signature(self, dec_places=3, separator='#'):
        r"""
        Returns the Vector signature in form #<d>.<ddd>#....
        """
        import re
        regexp = re.compile(r'[\s\[\]]+')

        with printoptions(precision=dec_places, suppress=True,
                          linewidth=1024, threshold=100000000):
            return separator.join(regexp.split(str(self)))

    def get_hash_signature(self, dec_places=3):
        r"""
        Returns  the sha256 hash of vector signature
        """
        import hashlib
        return hashlib.sha256(self.get_signature(dec_places)).hexdigest()

    def parallel_to_main(self):
        """Returns true if self is parallel to any main axis"""
        _l = self.leng
        _dots = numpy.abs(self.dot(Vector([1.0, 0.0, 0.0]))/_l)
        if numpy.allclose(_dots, 1.0):
            # all directions in self are parallel to X
            return True

        _dots += numpy.abs(self.dot(Vector([0.0, 1.0, 0.0]))/_l)
        if numpy.allclose(_dots, 1.0):
            # all directions in self are parallel to X or Y
            return True

        _dots += numpy.abs(self.dot(Vector([0.0, 0.0, 1.0]))/_l)
        if numpy.allclose(_dots, 1.0):
            # all directions in self are parallel to X or Y or Z
            return True

        return False

    def to_mayavi(self, **kw):
        r"""
        Returns a valid list of dictionaries for
        :func:`mayavirender.render <primfitter.formats.mayavirender.render>`.

        :class:`Vector` is rendered with
        :func:`mayavirender.points3d <primfitter.formats.mayavirender.points3d>`.

        By default the mesh is rendered without texture nor normals, but
        see keyword arguments.

        **Example**

        >>> from primfitter.geo import ufunc
        >>> from primfitter.formats import mayavirender
        >>> points, tris = ufunc.random_vectors(100, 10.0)

        Render the points, the fast way with default arguments:

        >>> pipe0 = mayavirender.render(points)

        This adds a *Point Cloud* to mayavi pipeline rendered in default
        *fgcolor*.

        Now tweak the pipeline by invoking the method :meth:`to_mayavi`
        with non default arguments:

        >>> lengs = points.leng
        >>> mlab_dict = points.to_mayavi(scalars=lengs, \
        ...     name="Points by Length", \
        ...     mode="sphere", \
        ...     scale_mode="scalar", \
        ...     scale_factor=0.1, \
        ...     sb_title="Lengths")
        >>> pipe1 = mayavirender.render(mlab_dict, figure="Points by Length", \
                size=(640, 480))

        **Keyword arguments**

        ``scalars`` *array shape(len(self))* ``default None``

            Array of scalar values associated to points,  will be colored
            using this value.

        ``name`` *str* ``default 'Point Cloud'``

            The name of the mlab_ *pipeline item* associated with the points.

        Any other keywords are added to the dictionary and passed to
        :func:`mayavirender.points3d <primfitter.formats.mayavirender.points3d>`
        without further checking. To see alist of default values use:

        >>> import primfitter.formats.mayavirender as mayavirender
        >>> print mayavirender.MAYAVI_CFG['common']
        >>> print mayavirender.MAYAVI_CFG['points3d']

        **Returns**

        ``[dict]``

            One valid dictionary for
            :func:`mayavirender.render <primfitter.formats.mayavirender.render>`.

        ..  seealso::

            The module :mod:`~primfitter.formats.mayavirender`.

        """
        _dict = {}
        _dict['args'] = [self]
        _dict['kw'] = {'name': "Point Cloud"}
        _dict['function'] = 'points3d'
        for _key in kw:
            if _key == 'scalars':
                if (isinstance(kw[_key], (list, tuple, numpy.ndarray))) \
                    and (len(kw[_key]) == len(self)):
                    _dict['args'].append(kw[_key])
                else:
                    print("Vector.to_mayavi: invalid scalars, ignored.",
                          file=sys.stderr)
            else:
                _dict['kw'][_key] = kw[_key]
        return [_dict]

    # Properties

    x = property(_get_x, _set_x)
    y = property(_get_y, _set_y)
    z = property(_get_z, _set_z)
    bbox = property(_get_bbox)
    bbox_center = property(_get_bbox_center)
    cyl = property(_get_cyl)
    cyl_deg = property(_get_cyl_deg)
    z_angle = property(_get_z_angle)
    z_angle_deg = property(_get_z_angle_deg)
    sp = property(_get_sp)
    sp_deg = property(_get_sp_deg)
    az = property(_get_az)
    az_deg = property(_get_az_deg)
    rad = property(_get_rad)
    el = property(_get_el)
    el_deg = property(_get_el_deg)
    leng = property(_get_leng)
    leng_sq = property(_get_leng_sq)
    unit = property(_get_unit)
    centroid = property(_get_centroid)
