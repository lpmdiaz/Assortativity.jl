# holds relevant node information for JSON export
struct JSON_node
   label::String
   group::Symbol
   group_index::Int
end

# holds relevant edge information for JSON export
struct JSON_edge
   source::String
   destination::String
   weight::Float64
end
