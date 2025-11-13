# required Tables.AbstractColumns object methods
# ## Tables interface
function Tables.getcolumn(d::AbstractStudyDesign, var::Symbol) # get a single variable
	if is_within(d, var)
		return getproperty(d.within, var)
	elseif is_between(d, var)
		btw_var = getproperty(d.between, var)
		return btw_var[d.uo.i]
	elseif is_covariate(d, var)
		return getproperty(d.covariates, var)
	else
		_err_not_in_design(var)
	end
end

Tables.istable(::Type{<:AbstractStudyDesign}) = true
Tables.columnaccess(::Type{<:AbstractStudyDesign}) = true
function Tables.columns(d::AbstractStudyDesign)::NamedTuple
	cov = isempty(d.covariates) ? (;) : Tables.columns(d.covariates)
	if d isa BetweenDesign
		return merge(_expand_between(d), cov)
	elseif d isa MixedDesign
		return merge(_expand_between(d), Tables.columns(d.within), cov)
	else
		# within design
		return merge((; d.uo.name => d.uo.values), Tables.columns(d.within), cov)
	end
end
Tables.getcolumn(d::AbstractStudyDesign, ::Type{T}, col::Int, var::Symbol) where {T} = getcolumn(d, var)
Tables.getcolumn(d::AbstractStudyDesign, i::Int) = getcolumn(d, names(d)[i])
Tables.columnnames(d::AbstractStudyDesign) = names(d)

DataAPI.nrow(d::AbstractStudyDesign) = length(d.uo.i)
DataAPI.ncol(d::AbstractStudyDesign) = length(names(d))


@inline function _expand_between(d::AbstractStudyDesign)::NamedTuple
	if d isa WithinDesign
		return (;)
	else
		# expand between design to full length with unit of observations
		return select_rows(columntable(d.between), d.uo.i)
	end
end
