var documenterSearchIndex = {"docs":
[{"location":"api/#API-Reference","page":"API Reference","title":"API Reference","text":"","category":"section"},{"location":"api/#VacuumFields","page":"API Reference","title":"VacuumFields","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"AbstractCircuit\nAbstractCoil\nAbstractControlPoint\nDistributedCoil\nFluxControlPoint\nFluxControlPoints\nGS_IMAS_pf_active__coil\nIMAS_pf_active__coils\nIsoControlPoint\nIsoControlPoints\nMultiCoil\nParallelogramCoil\nPointCoil\nQuadCoil\nSaddleControlPoint\nSaddleControlPoints\nSeriesCircuit\nboundary_control_points\nd2M_dZ2\ndM_dZ\nencircling_coils\nfind_coil_currents!\nfixed2free\nmutual\nnormalized_growth_rate\noptimal_λ_regularize\nstability_margin\nupdate_coil_currents!","category":"page"},{"location":"api/#VacuumFields.AbstractCoil","page":"API Reference","title":"VacuumFields.AbstractCoil","text":"AbstractCoil\n\nAbstract coil type\n\n\n\n\n\n","category":"type"},{"location":"api/#VacuumFields.AbstractControlPoint","page":"API Reference","title":"VacuumFields.AbstractControlPoint","text":"AbstractControlPoint{T<:Real}\n\nAbstract control point - for coil current least-squres fitting\n\n\n\n\n\n","category":"type"},{"location":"api/#VacuumFields.DistributedCoil","page":"API Reference","title":"VacuumFields.DistributedCoil","text":"DistributedCoil{T1<:Real,T2<:Real,T3<:Real,T4<:Real} <: AbstractSingleCoil{T1,T2,T3,T4}\n\nCoil consisting of distributed filaments\n\n\n\n\n\n","category":"type"},{"location":"api/#VacuumFields.FluxControlPoint","page":"API Reference","title":"VacuumFields.FluxControlPoint","text":"FluxControlPoint(R::Real, Z::Real, target::Real, weight::Real=1.0)\n\nReturn a control point for a target flux value at point (R, Z), with an optional weight\n\n\n\n\n\n","category":"type"},{"location":"api/#VacuumFields.FluxControlPoints","page":"API Reference","title":"VacuumFields.FluxControlPoints","text":"FluxControlPoints(Rs::AbstractVector{<:Real}, Zs::AbstractVector{<:Real}, ψtarget::Real)\n\nReturn a Vector of FluxControlPoint at each Rs and Zs point, each with the same ψtarget flux\n\n\n\n\n\nFluxControlPoints(Rs::AbstractVector{<:Real}, Zs::AbstractVector{<:Real}, ψtargets::AbstractVector{<:Real})\n\nReturn a Vector of FluxControlPoint at each Rs and Zs point with ψtargets flux\n\n\n\n\n\n","category":"function"},{"location":"api/#VacuumFields.IsoControlPoint","page":"API Reference","title":"VacuumFields.IsoControlPoint","text":"IsoControlPoint(R::Real, Z::Real,weight::Real=1.0)\n\nReturn a control point for equal flux between points (R1, Z1) and (R2, Z2), with an optional weight\n\n\n\n\n\n","category":"type"},{"location":"api/#VacuumFields.IsoControlPoints","page":"API Reference","title":"VacuumFields.IsoControlPoints","text":"IsoControlPoints(Rs::AbstractVector{<:Real}, Zs::AbstractVector{<:Real})\n\nReturn a Vector of IsoControlPoints between each pair of Rs and Zs points\n\n\n\n\n\n","category":"function"},{"location":"api/#VacuumFields.MultiCoil","page":"API Reference","title":"VacuumFields.MultiCoil","text":"MultiCoil{SC<:AbstractSingleCoil} <: AbstractCoil\n\nA coil consisting of multiple coils linke in series\n\n\n\n\n\n","category":"type"},{"location":"api/#VacuumFields.ParallelogramCoil","page":"API Reference","title":"VacuumFields.ParallelogramCoil","text":"ParallelogramCoil{T1, T2, T3, T4} <:  AbstractSingleCoil{T1, T2, T3, T4}\n\nParallelogram coil with the R, Z, ΔR, ΔZ, θ1, θ2 formalism (as used by EFIT, for example). Here θ1 and θ2 are the shear angles along the x- and y-axes, respectively, in degrees.\n\n\n\n\n\n","category":"type"},{"location":"api/#VacuumFields.PointCoil","page":"API Reference","title":"VacuumFields.PointCoil","text":"PointCoil{T1, T2, T3, T4} <:  AbstractSingleCoil{T1, T2, T3, T4}\n\nPoint filament coil at scalar (R, Z) location\n\n\n\n\n\n","category":"type"},{"location":"api/#VacuumFields.QuadCoil","page":"API Reference","title":"VacuumFields.QuadCoil","text":"QuadCoil{T1, T2, T3, T4} <:  AbstractSingleCoil{T1, T2, T3, T4}\n\nQuadrilateral coil with counter-clockwise corners (starting from lower left) at R and Z\n\n\n\n\n\n","category":"type"},{"location":"api/#VacuumFields.SaddleControlPoint","page":"API Reference","title":"VacuumFields.SaddleControlPoint","text":"SaddleControlPoint(R::Real, Z::Real,weight::Real=1.0)\n\nReturn a control point for a saddle (i.e., dψ/dR = dψ/dZ = 0.0) at point (R, Z), with an optional weight\n\n\n\n\n\n","category":"type"},{"location":"api/#VacuumFields.SaddleControlPoints","page":"API Reference","title":"VacuumFields.SaddleControlPoints","text":"SaddleControlPoints(Rs::AbstractVector{<:Real}, Zs::AbstractVector{<:Real})\n\nReturn a Vector of SaddleControlPoint at each Rs and Zs point\n\n\n\n\n\n","category":"function"},{"location":"api/#VacuumFields.boundary_control_points","page":"API Reference","title":"VacuumFields.boundary_control_points","text":"boundary_control_points(EQfixed::MXHEquilibrium.AbstractEquilibrium, fraction_inside::Float64=0.999, ψbound::Real=0.0; Npts::Integer=99)\n\nReturn a Vector of FluxControlPoint, each with target ψbound, at Npts equally distributed fraction_inside percent inside the the boundary of EQfixed\n\n\n\n\n\n","category":"function"},{"location":"api/#VacuumFields.d2M_dZ2","page":"API Reference","title":"VacuumFields.d2M_dZ2","text":"d2M_dZ2(EQ::MXHEquilibrium.AbstractEquilibrium, C::Union{AbstractCoil, IMAScoil}, δZ::Real=0.0;\n      COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(EQ), kwargs...)\n\nCompute the second Z derivative of the mutual inductance between an equilibrium and a coil,     where the equilibrium is shifted vertically by δZ\n\n\n\n\n\nd2M_dZ2(image::Image, C::Union{AbstractCoil, IMAScoil}, Ip::Real, δZ::Real=0.0;\n      COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(EQ), kwargs...)\n\nCompute the second Z derivative of the mutual inductance between an equilibrium's image current and a coil,     where the equilibrium is shifted vertically by δZ\n\n\n\n\n\n","category":"function"},{"location":"api/#VacuumFields.dM_dZ","page":"API Reference","title":"VacuumFields.dM_dZ","text":"dM_dZ(EQ::MXHEquilibrium.AbstractEquilibrium, C::Union{AbstractCoil, IMAScoil}, δZ::Real=0.0;\n      COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(EQ), kwargs...)\n\nCompute the Z derivative of the mutual inductance between an equilibrium and a coil,     where the equilibrium is shifted vertically by δZ\n\n\n\n\n\ndM_dZ(image::Image, C::Union{AbstractCoil, IMAScoil}, Ip::Real, δZ::Real=0.0;\n      COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(EQ), kwargs...)\n\nCompute the Z derivative of the mutual inductance between an equilibrium's image current and a coil,     where the equilibrium is shifted vertically by δZ\n\n\n\n\n\n","category":"function"},{"location":"api/#VacuumFields.encircling_coils","page":"API Reference","title":"VacuumFields.encircling_coils","text":"encircling_coils(EQfixed::MXHEquilibrium.AbstractEquilibrium, n_coils::Int)\n\nReturn a Vector of n_coils PointCoils distributed outside of EQfixed's boundary\n\n\n\n\n\nencircling_coils(bnd_r::AbstractVector{T}, bnd_z::AbstractVector{T}, n_coils::Int) where {T<:Real}\n\nReturn a Vector of n_coils PointCoils distributed outside of closed boundary defined by bnd_r and bnd_z\n\n\n\n\n\n","category":"function"},{"location":"api/#VacuumFields.find_coil_currents!","page":"API Reference","title":"VacuumFields.find_coil_currents!","text":"find_coil_currents!(\n    coils::AbstractVector{<:AbstractCoil},\n    EQ::MXHEquilibrium.AbstractEquilibrium;\n    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],\n    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],\n    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],\n    ψbound::Real=0.0,\n    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],\n    λ_regularize::Float64=0.0,\n    Sb::MXHEquilibrium.Boundary=plasma_boundary_psi_w_fallback(EQ)[1])\n\nFind the currents for coils that best match (least-squares) the control points provided by flux_cps, saddle_cps, and iso_cps\n\nAssumes flux from  plasma current given by equilibrium EQ with a ψbound flux at the boundary Sb\n\nOptionally assumes flux from additional fixed_coils, whose currents will not change\n\nλ_regularize provides regularization in the least-squares fitting\n\n\n\n\n\nfind_coil_currents!(\n    coils::AbstractVector{<:AbstractCoil},\n    EQ::MXHEquilibrium.AbstractEquilibrium,\n    image::Image;\n    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],\n    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],\n    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],\n    ψbound::Real=0.0,\n    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],\n    λ_regularize::Float64=0.0,\n    Sb::MXHEquilibrium.Boundary=plasma_boundary_psi_w_fallback(EQ)[1])\n\nFind the currents for coils that best match (leas-squares) the control points provided by flux_cps, saddle_cps, and iso_cps\n\nAssumes flux from  plasma current given by equilibrium EQ with image currents image and a ψbound flux at the boundary Sb\n\nOptionally assumes flux from additional fixed_coils, whose currents will not change\n\nλ_regularize provides regularization in the least-squares fitting\n\n\n\n\n\nfind_coil_currents!(\n    coils::AbstractVector{<:AbstractCoil},\n    EQ::MXHEquilibrium.AbstractEquilibrium,\n    image::Image;\n    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],\n    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],\n    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],\n    ψbound::Real=0.0,\n    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],\n    λ_regularize::Float64=0.0,\n    Sb::MXHEquilibrium.Boundary=plasma_boundary_psi_w_fallback(EQ)[1])\n\nFind the currents for coils that best match (leas-squares) the control points provided by flux_cps, saddle_cps, and iso_cps\n\nAssumes flux from  plasma current given by equilibrium EQ with image currents image and a ψbound flux at the boundary Sb\n\nOptionally assumes flux from additional fixed_coils, whose currents will not change\n\nλ_regularize provides regularization in the least-squares fitting\n\n\n\n\n\nfind_coil_currents!(\n    coils::AbstractVector{<:AbstractCoil},\n    EQ::Nothing;\n    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],\n    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],\n    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],\n    ψbound::Real=0.0,\n    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],\n    λ_regularize::Float64=0.0,\n    cocos=MXHEquilibrium.cocos(11),\n    Sb=nothing)\n\nFind the currents for coils that best match (leas-squares) the control points provided by flux_cps, saddle_cps, and iso_cps\n\nVacuume case: assumes no equilibrium plasma current\n\nOptionally assumes flux from additional fixed_coils, whose currents will not change\n\nλ_regularize provides regularization in the least-squares fitting\n\n\n\n\n\nfind_coil_currents!(\n    coils::AbstractVector{<:AbstractCoil},\n    EQ::Nothing, # VACUUM case\n    image::Nothing;\n    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],\n    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],\n    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],\n    ψbound::Real=0.0,\n    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],\n    λ_regularize::Float64=0.0,\n    cocos=MXHEquilibrium.cocos(11),\n    Sb=nothing)\n\nFind the currents for coils that best match (leas-squares) the control points provided by flux_cps, saddle_cps, and iso_cps\n\nVacuume case: assumes no equilibrium plasma current\n\nOptionally assumes flux from additional fixed_coils, whose currents will not change\n\nλ_regularize provides regularization in the least-squares fitting\n\n\n\n\n\nfind_coil_currents!(\n    coils::AbstractVector{<:AbstractCoil};\n    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],\n    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],\n    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],\n    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],\n    λ_regularize::Float64=0.0,\n    cocos=MXHEquilibrium.cocos(11))\n\nFind the currents for coils that best match (leas-squares) the control points provided by flux_cps, saddle_cps, and iso_cps\n\nVacuume case: assumes no equilibrium plasma current\n\nOptionally assumes flux from additional fixed_coils, whose currents will not change\n\nλ_regularize provides regularization in the least-squares fitting\n\n\n\n\n\n","category":"function"},{"location":"api/#VacuumFields.fixed2free","page":"API Reference","title":"VacuumFields.fixed2free","text":"fixed2free(\n    EQfixed::MXHEquilibrium.AbstractEquilibrium,\n    coils::AbstractVector{<:AbstractCoil},\n    R::AbstractVector{T},\n    Z::AbstractVector{T};\n    ψbound::Real=0.0,\n    kwargs...) where {T<:Real}\n\nConvert the flux of a fixed-boundary equilibrium EQfixed to a free-boundary representation on an (R,Z) grid, using the flux from coils with currents satisfying given control points\n\n\n\n\n\nfixed2free(\n    EQfixed::MXHEquilibrium.AbstractEquilibrium,\n    image::Image,\n    coils::AbstractVector{<:AbstractCoil},\n    R::AbstractVector{T},\n    Z::AbstractVector{T};\n    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],\n    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],\n    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],\n    ψbound::Real=0.0,\n    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],\n    λ_regularize::Real=0.0) where {T<:Real}\n\nConvert the flux of a fixed-boundary equilibrium EQfixed to a free-boundary representation on an (R,Z) grid, subtracting out the flux from image currents image and adding in the flux from coils with currents that best satisfy the control points and fixed_coils with predefined coil currents\n\n\n\n\n\n","category":"function"},{"location":"api/#VacuumFields.mutual","page":"API Reference","title":"VacuumFields.mutual","text":"mutual(C1::Union{AbstractCoil, IMASelement}, C2::PointCoil)\n\nCompute the mutual inductance between an arbitrary coil or IMAS.pf_active__coil___element and a PointCoil\n\n\n\n\n\nmutual(C1::Union{AbstractCoil, IMASelement}, C2::DistributedCoil)\n\nCompute the mutual inductance between an arbitrary coil or IMAS.pf_active__coil___element and a DistributedCoil\n\n\n\n\n\nmutual(C1::Union{AbstractCoil, IMASelement}, C2::Union{ParallelogramCoil, QuadCoil, IMASelement}; xorder::Int=3, yorder::Int=3)\n\nCompute the mutual inductance between an arbitrary coil or IMAS.pf_active__coil___element and     a ParallelogramCoil, QuadCoil, or IMAS.pf_active__coil___element\n\nxorder and yorder give the order of Gauss-Legendre quadrature for integration over the coil area\n\n\n\n\n\nmutual(C1::Union{AbstractCoil, IMASelement}, C2::IMAScoil; xorder::Int=3, yorder::Int=3)\n\nCompute the mutual inductance between an arbitrary coil or IMAS.pf_active__coil___element and a IMAS.pf_active__coil\n\nxorder and yorder give the order of Gauss-Legendre quadrature for integration over the coil area\n\n\n\n\n\nmutual(C1::IMAScoil, C2::Union{AbstractCoil, IMAScoil, IMASelement}; xorder::Int=3, yorder::Int=3)\n\nCompute the mutual inductance between an IMAS.pf_active__coil and an arbitrary coil, IMAS.pf_active__coil___element, or a IMAS.pf_active__coil\n\nxorder and yorder give the order of Gauss-Legendre quadrature for integration over the coil area\n\n\n\n\n\nmutual(EQ::MXHEquilibrium.AbstractEquilibrium, C::Union{AbstractCoil, IMAScoil}, δZ::Real=0.0;\n    COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(EQ), kwargs...)\n\nCompute the mutual inductance between an equilibrium and a coil,     where the equilibrium is shifted vertically by δZ\n\n\n\n\n\nmutual(image::Image, C::Union{AbstractCoil, IMAScoil}, Ip::Real, δZ::Real=0.0;\n       COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11), kwargs...)\n\nCompute the mutual inductance between an equilibrium's image current and a coil,     where the equilibrium is shifted vertically by δZ\n\n\n\n\n\n","category":"function"},{"location":"api/#VacuumFields.normalized_growth_rate","page":"API Reference","title":"VacuumFields.normalized_growth_rate","text":"normalized_growth_rate(EQ::MXHEquilibrium.AbstractEquilibrium, coils::Vector{<:Union{AbstractCoil, IMAScoil}}; kwargs...)\n\nCompute the vertical growth rate (γ) and effective vertical time constant (τ, weighted L/R time) for a given equilibrium and coils\n\nReturn (γ, τ, γ * τ), where γ * τ < 10 for stability or controllability\n\nThis is the massless approximation and only use the passive conductors for computing τ (per advice from Olofsson)\n\n\n\n\n\nnormalized_growth_rate(image::Image, coils::Vector{<:Union{AbstractCoil, IMAScoil}}, Ip::Real; order::Int=default_order)\n\nCompute the vertical growth rate (γ) and effective vertical time constant (τ, weighted L/R time), for a given equilibrium's image & plasma current and coils\n\nReturn (γ, τ, γ * τ), where γ * τ < 10 for stability or controllability\n\nThis is the massless approximation and only use the passive conductors for computing τ (per advice from Olofsson)\n\n\n\n\n\n","category":"function"},{"location":"api/#VacuumFields.optimal_λ_regularize","page":"API Reference","title":"VacuumFields.optimal_λ_regularize","text":"optimal_λ_regularize(\n    coils::AbstractVector{<:AbstractCoil},\n    EQ::MXHEquilibrium.AbstractEquilibrium,\n    image::Image;\n    flux_cps::Vector{<:FluxControlPoint}=FluxControlPoint{Float64}[],\n    saddle_cps::Vector{<:SaddleControlPoint}=SaddleControlPoint{Float64}[],\n    iso_cps::Vector{<:IsoControlPoint}=IsoControlPoint{Float64}[],\n    ψbound::Real=0.0,\n    fixed_coils::AbstractVector{<:AbstractCoil}=PointCoil{Float64,Float64}[],\n    min_exp::Integer=-20, max_exp::Integer=-10,\n    Sb::MXHEquilibrium.Boundary=plasma_boundary_psi_w_fallback(EQ)[1])\n\nFind the optimizal λ_regularize to be used within find_coil_currents! that minimizes the cost\n\n\n\n\n\n","category":"function"},{"location":"api/#VacuumFields.stability_margin","page":"API Reference","title":"VacuumFields.stability_margin","text":"stability_margin(EQ::MXHEquilibrium.AbstractEquilibrium, coils::Vector{<:Union{AbstractCoil, IMAScoil}}; kwargs...)\n\nCompute the m_s inductive stability margin for a given equilibrium and coils.     Should be greater than 0.15 for vertical stability\n\nFirst introduced in A. Portone, Nucl. Fusion 45 (2005) 926–932. https://doi.org/10.1088/0029-5515/45/8/021\n\n\n\n\n\nstability_margin(image::Image, coils::Vector{<:Union{AbstractCoil, IMAScoil}}, Ip::Real; order::Int=default_order)\n\nCompute the m_s inductive stability margin for a given equilibrium's image & plasma current and coils.     Should be greater than 0.15 for vertical stability\n\nFirst introduced in A. Portone, Nucl. Fusion 45 (2005) 926–932. https://doi.org/10.1088/0029-5515/45/8/021\n\n\n\n\n\n","category":"function"},{"location":"license/","page":"License","title":"License","text":"                             Apache License\n                       Version 2.0, January 2004\n                    http://www.apache.org/licenses/","category":"page"},{"location":"license/","page":"License","title":"License","text":"TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION","category":"page"},{"location":"license/","page":"License","title":"License","text":"Definitions.\n\"License\" shall mean the terms and conditions for use, reproduction, and distribution as defined by Sections 1 through 9 of this document.\n\"Licensor\" shall mean the copyright owner or entity authorized by the copyright owner that is granting the License.\n\"Legal Entity\" shall mean the union of the acting entity and all other entities that control, are controlled by, or are under common control with that entity. For the purposes of this definition, \"control\" means (i) the power, direct or indirect, to cause the direction or management of such entity, whether by contract or otherwise, or (ii) ownership of fifty percent (50%) or more of the outstanding shares, or (iii) beneficial ownership of such entity.\n\"You\" (or \"Your\") shall mean an individual or Legal Entity exercising permissions granted by this License.\n\"Source\" form shall mean the preferred form for making modifications, including but not limited to software source code, documentation source, and configuration files.\n\"Object\" form shall mean any form resulting from mechanical transformation or translation of a Source form, including but not limited to compiled object code, generated documentation, and conversions to other media types.\n\"Work\" shall mean the work of authorship, whether in Source or Object form, made available under the License, as indicated by a copyright notice that is included in or attached to the work (an example is provided in the Appendix below).\n\"Derivative Works\" shall mean any work, whether in Source or Object form, that is based on (or derived from) the Work and for which the editorial revisions, annotations, elaborations, or other modifications represent, as a whole, an original work of authorship. For the purposes of this License, Derivative Works shall not include works that remain separable from, or merely link (or bind by name) to the interfaces of, the Work and Derivative Works thereof.\n\"Contribution\" shall mean any work of authorship, including the original version of the Work and any modifications or additions to that Work or Derivative Works thereof, that is intentionally submitted to Licensor for inclusion in the Work by the copyright owner or by an individual or Legal Entity authorized to submit on behalf of the copyright owner. For the purposes of this definition, \"submitted\" means any form of electronic, verbal, or written communication sent to the Licensor or its representatives, including but not limited to communication on electronic mailing lists, source code control systems, and issue tracking systems that are managed by, or on behalf of, the Licensor for the purpose of discussing and improving the Work, but excluding communication that is conspicuously marked or otherwise designated in writing by the copyright owner as \"Not a Contribution.\"\n\"Contributor\" shall mean Licensor and any individual or Legal Entity on behalf of whom a Contribution has been received by Licensor and subsequently incorporated within the Work.\nGrant of Copyright License. Subject to the terms and conditions of this License, each Contributor hereby grants to You a perpetual, worldwide, non-exclusive, no-charge, royalty-free, irrevocable copyright license to reproduce, prepare Derivative Works of, publicly display, publicly perform, sublicense, and distribute the Work and such Derivative Works in Source or Object form.\nGrant of Patent License. Subject to the terms and conditions of this License, each Contributor hereby grants to You a perpetual, worldwide, non-exclusive, no-charge, royalty-free, irrevocable (except as stated in this section) patent license to make, have made, use, offer to sell, sell, import, and otherwise transfer the Work, where such license applies only to those patent claims licensable by such Contributor that are necessarily infringed by their Contribution(s) alone or by combination of their Contribution(s) with the Work to which such Contribution(s) was submitted. If You institute patent litigation against any entity (including a cross-claim or counterclaim in a lawsuit) alleging that the Work or a Contribution incorporated within the Work constitutes direct or contributory patent infringement, then any patent licenses granted to You under this License for that Work shall terminate as of the date such litigation is filed.\nRedistribution. You may reproduce and distribute copies of the Work or Derivative Works thereof in any medium, with or without modifications, and in Source or Object form, provided that You meet the following conditions:\n(a) You must give any other recipients of the Work or     Derivative Works a copy of this License; and\n(b) You must cause any modified files to carry prominent notices     stating that You changed the files; and\n(c) You must retain, in the Source form of any Derivative Works     that You distribute, all copyright, patent, trademark, and     attribution notices from the Source form of the Work,     excluding those notices that do not pertain to any part of     the Derivative Works; and\n(d) If the Work includes a \"NOTICE\" text file as part of its     distribution, then any Derivative Works that You distribute must     include a readable copy of the attribution notices contained     within such NOTICE file, excluding those notices that do not     pertain to any part of the Derivative Works, in at least one     of the following places: within a NOTICE text file distributed     as part of the Derivative Works; within the Source form or     documentation, if provided along with the Derivative Works; or,     within a display generated by the Derivative Works, if and     wherever such third-party notices normally appear. The contents     of the NOTICE file are for informational purposes only and     do not modify the License. You may add Your own attribution     notices within Derivative Works that You distribute, alongside     or as an addendum to the NOTICE text from the Work, provided     that such additional attribution notices cannot be construed     as modifying the License.\nYou may add Your own copyright statement to Your modifications and may provide additional or different license terms and conditions for use, reproduction, or distribution of Your modifications, or for any such Derivative Works as a whole, provided Your use, reproduction, and distribution of the Work otherwise complies with the conditions stated in this License.\nSubmission of Contributions. Unless You explicitly state otherwise, any Contribution intentionally submitted for inclusion in the Work by You to the Licensor shall be under the terms and conditions of this License, without any additional terms or conditions. Notwithstanding the above, nothing herein shall supersede or modify the terms of any separate license agreement you may have executed with Licensor regarding such Contributions.\nTrademarks. This License does not grant permission to use the trade names, trademarks, service marks, or product names of the Licensor, except as required for reasonable and customary use in describing the origin of the Work and reproducing the content of the NOTICE file.\nDisclaimer of Warranty. Unless required by applicable law or agreed to in writing, Licensor provides the Work (and each Contributor provides its Contributions) on an \"AS IS\" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied, including, without limitation, any warranties or conditions of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A PARTICULAR PURPOSE. You are solely responsible for determining the appropriateness of using or redistributing the Work and assume any risks associated with Your exercise of permissions under this License.\nLimitation of Liability. In no event and under no legal theory, whether in tort (including negligence), contract, or otherwise, unless required by applicable law (such as deliberate and grossly negligent acts) or agreed to in writing, shall any Contributor be liable to You for damages, including any direct, indirect, special, incidental, or consequential damages of any character arising as a result of this License or out of the use or inability to use the Work (including but not limited to damages for loss of goodwill, work stoppage, computer failure or malfunction, or any and all other commercial damages or losses), even if such Contributor has been advised of the possibility of such damages.\nAccepting Warranty or Additional Liability. While redistributing the Work or Derivative Works thereof, You may choose to offer, and charge a fee for, acceptance of support, warranty, indemnity, or other liability obligations and/or rights consistent with this License. However, in accepting such obligations, You may act only on Your own behalf and on Your sole responsibility, not on behalf of any other Contributor, and only if You agree to indemnify, defend, and hold each Contributor harmless for any liability incurred by, or claims asserted against, such Contributor by reason of your accepting any such warranty or additional liability.","category":"page"},{"location":"license/","page":"License","title":"License","text":"END OF TERMS AND CONDITIONS","category":"page"},{"location":"license/","page":"License","title":"License","text":"APPENDIX: How to apply the Apache License to your work.","category":"page"},{"location":"license/","page":"License","title":"License","text":"  To apply the Apache License to your work, attach the following\n  boilerplate notice, with the fields enclosed by brackets \"[]\"\n  replaced with your own identifying information. (Don't include\n  the brackets!)  The text should be enclosed in the appropriate\n  comment syntax for the file format. We also recommend that a\n  file or class name and description of purpose be included on the\n  same \"printed page\" as the copyright notice for easier\n  identification within third-party archives.","category":"page"},{"location":"license/","page":"License","title":"License","text":"Copyright 2024 General Atomics","category":"page"},{"location":"license/","page":"License","title":"License","text":"Licensed under the Apache License, Version 2.0 (the \"License\");    you may not use this file except in compliance with the License.    You may obtain a copy of the License at","category":"page"},{"location":"license/","page":"License","title":"License","text":"   http://www.apache.org/licenses/LICENSE-2.0","category":"page"},{"location":"license/","page":"License","title":"License","text":"Unless required by applicable law or agreed to in writing, software    distributed under the License is distributed on an \"AS IS\" BASIS,    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    See the License for the specific language governing permissions and    limitations under the License.","category":"page"},{"location":"#VacuumFields.jl","page":"VacuumFields.jl","title":"VacuumFields.jl","text":"","category":"section"},{"location":"","page":"VacuumFields.jl","title":"VacuumFields.jl","text":"Closed-boundary solvers do not provide a poloidal field solution in the vacuum region, outside of the last closed flux surface. VacuumFields.jl uses the Green’s function method to match the plasma current contribution with the contributions of external poloidal field coils at the last closed flux surface to determine the homogeneous solution. The total solution is then extended into the vacuum region to get a realistic vacuum solution.","category":"page"},{"location":"#Online-documentation","page":"VacuumFields.jl","title":"Online documentation","text":"","category":"section"},{"location":"","page":"VacuumFields.jl","title":"VacuumFields.jl","text":"For more details, see the online documentation.","category":"page"},{"location":"","page":"VacuumFields.jl","title":"VacuumFields.jl","text":"(Image: Docs)","category":"page"},{"location":"notice/#VacuumFields.jl-Notice","page":"Notice","title":"VacuumFields.jl Notice","text":"","category":"section"},{"location":"notice/","page":"Notice","title":"Notice","text":"The purpose of this NOTICE file is to provide legal notices and acknowledgments that must be displayed to users in any derivative works or distributions. This file does not alter the terms of the Apache 2.0 license that governs the use and distribution of the VacuumFields.jl package.","category":"page"},{"location":"notice/#Development-Attribution","page":"Notice","title":"Development Attribution","text":"","category":"section"},{"location":"notice/","page":"Notice","title":"Notice","text":"VacuumFields.jl was originally developed under the FUSE project by the Magnetic Fusion Energy group at General Atomics.","category":"page"},{"location":"notice/#Citation","page":"Notice","title":"Citation","text":"","category":"section"},{"location":"notice/","page":"Notice","title":"Notice","text":"If this software contributes to an academic publication, please cite it as follows:","category":"page"},{"location":"notice/","page":"Notice","title":"Notice","text":"@article{meneghini2024fuse,\nauthor = {Meneghini, O. and Slendebroek, T. and Lyons, B.C. and McLaughlin, K. and McClenaghan, J. and Stagner, L. and Harvey, J. and Neiser, T.F. and Ghiozzi, A. and Dose, G. and Guterl, J. and Zalzali, A. and Cote, T. and Shi, N. and Weisberg, D. and Smith, S.P. and Grierson, B.A. and Candy, J.},\ndoi = {10.48550/arXiv.2409.05894},\njournal = {arXiv},\ntitle = {{FUSE (Fusion Synthesis Engine): A Next Generation Framework for Integrated Design of Fusion Pilot Plants}},\nyear = {2024}\n}","category":"page"},{"location":"notice/#Trademark-Notice","page":"Notice","title":"Trademark Notice","text":"","category":"section"},{"location":"notice/","page":"Notice","title":"Notice","text":"The names \"General Atomics\", and any associated logos or images, are trademarks of General Atomics. Use of these trademarks without prior written consent from General Atomics is strictly prohibited. Users cannot imply endorsement by General Atomics or contributors to the project simply because the project is part of their work.","category":"page"},{"location":"notice/#Copyright","page":"Notice","title":"Copyright","text":"","category":"section"},{"location":"notice/","page":"Notice","title":"Notice","text":"Copyright (c) 2024 General Atomics","category":"page"},{"location":"notice/#Version","page":"Notice","title":"Version","text":"","category":"section"},{"location":"notice/","page":"Notice","title":"Notice","text":"Version: v2.1","category":"page"}]
}
