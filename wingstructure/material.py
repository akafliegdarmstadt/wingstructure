from abc import ABC, abstractmethod

class _AbstractMaterial(ABC):
    def __init__(self, ρ):
        self.ρ = ρ
    

class IsotropicMaterial(_AbstractMaterial):
    def __init__(self, ρ, E, G):
        super().__init__(ρ)
        self.E = E
        self.G = G
    
    