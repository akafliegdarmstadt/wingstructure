"""This module contains the Masspoint class for center of gravity calculations"""


class Masspoint(object):
    """Masspoints can be added to get the sum of masses and the common center of gravity"""

    def __init__(self, mass: float, point: tuple):
        """
        constructor for Masspoint
        :param mass: masspoints mass
        :param point: center of gravity of masspoint
        """
        self.mass = mass
        self.point = point

    def __add__(self, othermasspoint):
        """
        Calculates sum of masses and common center of gravity
        :param othermasspoint: another instance of Masspoint
        :return: new instance of Masspoint, with sum of mass and common cg
        """
        mass = self.mass + othermasspoint.mass
        cg = tuple(
            (self.mass * cg1 + othermasspoint.mass * cg2) \
            / mass for cg1, cg2 in zip(self.point, othermasspoint.point))

        return Masspoint(mass, cg)

    def __radd__(self, othermasspoint):
        """
        Calculates sum of masses and common center of gravity
        Is needed, when sum of Masspoints is calculated
        :param othermasspoint: another instance of Masspoint
        :return: new instance of Masspoint, with sum of mass and common cg
        """
        if not isinstance(othermasspoint, Masspoint):
            return self

    def __str__(self):
        return 'm={} at {}'.format(self.mass, self.point)

    def __repr__(self):
        return self.__str__()