"""
    AbstractDomain{S<:Union{AbstractFloat,arb},T<:Union{AbstractFloat,arb,Rational,fmpq}}

Abstract supertype for all domains. Most values of the domain are
given in the type `S`, the exception is values corresponding to angles
which are given in type `T`. If `T` is `Rational` or `fmpq` then the
angle values do not represent the angles directly but after a division
by `π`, this allows representing angles which are rational multiples
of `π` exactly.
"""
abstract type AbstractDomain{
    S<:Union{AbstractFloat,arb},
    T<:Union{AbstractFloat,arb,Rational,fmpq},
} end
