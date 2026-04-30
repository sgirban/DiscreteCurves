module AbstractTypes

export AbstractDiscreteCurve, TopologyTrait, ClosedTopology, OpenTopology, StorageTrait, DenseStorage, SparseStorage, CurvatureTrait, HasCurvatureData, NoCurvatureData, OrientationTrait, HasOrientationData, NoOrientationData
export topology_trait, is_closed_topology, is_open_topology, storage_trait, curvature_trait, orientation_trait

abstract type AbstractDiscreteCurve{N, T <: AbstractFloat} end

abstract type TopologyTrait end 
struct ClosedTopology <: TopologyTrait end
struct OpenTopology <: TopologyTrait end

topology_trait(::AbstractDiscreteCurve) = OpenTopology()
topology_trait(::Type{<:AbstractDiscreteCurve}) = OpenTopology()

@inline is_closed_topology(curve) = topology_trait(curve) isa ClosedTopology
@inline is_open_topology(curve) = topology_trait(curve) isa OpenTopology

abstract type StorageTrait end
abstract type DenseStorage <: StorageTrait end
abstract type SparseStorage <: StorageTrait end

storage_trait(::AbstractDiscreteCurve) = DenseStorage()

abstract type CurvatureTrait end
struct HasCurvatureData <: CurvatureTrait end
struct NoCurvatureData <: CurvatureTrait end

curvature_trait(::AbstractDiscreteCurve) = NoCurvatureData()

abstract type OrientationTrait end
struct HasOrientationData <: OrientationTrait end
struct NoOrientationData <: OrientationTrait end

orientation_trait(::AbstractDiscreteCurve) = NoOrientationData()

end