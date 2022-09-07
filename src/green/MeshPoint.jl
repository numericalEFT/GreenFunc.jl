"""
Cartisian product of 1 dimensional meshes 
"""

"""

#Members:
"""
struct MeshPoint{T}
    linear_index::Int
    index::Array{::Int,1}
    value::T
    weight::T
    function MeshPoint{T}(linear_index::Int, index::Int, value::T, weight::T
    ) 
        mnew = new{T}(linear_index,index,value,weight
        )
        return mnew
    end
end


"""
def __getitem__(self, *args):
        return self.value.__getitem__(*args)
"""
function Base.getindex(obj::MeshPoint, key...)

end


function Base.:-(obj::MeshPoint)

end

function Base.:+(obj::MeshPoint)

end

function Base.:+(obj1::MeshPoint, obj2::MeshPoint)

end

function Base.:-(obj1::MeshPoint, obj2::MeshPoint)

end

function Base.:*(obj1::MeshPoint, obj2::MeshPoint)

end

function Base.:/(obj1::MeshPoint, obj2::MeshPoint)

end


function Base.show(object::MeshPoint)

end
