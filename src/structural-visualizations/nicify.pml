set cartoon_side_chain_helper, on
#set float_labels, on
set label_size, 25
set label_font_id, 7
set ray_trace_mode, 0
set label_color,black
set label_padding, [0, 0, 0]
set label_bg_transparency, 0.2
set label_bg_color, white
center
disable
bg_color white
remove hydrogens
set label_position, (0, 0, 20)
set cartoon_transparency, 0.7
set cartoon_fancy_helices

# From Preben
# reinitialize
#set ray_opaque_background, off
#cmd.space('cmyk')
# viewport 400,2000
#util.cbaw
#bg white
#set light_count,8
#set spec_count,1
#set shininess, 10
#set specular, 0.25
#set ambient,0
#set direct,1
#set reflect,1.5
#set ray_shadow_decay_factor, 0.1
#set ray_shadow_decay_range, 2
#unset depth_cue
#set ray_trace_mode,3
#set sphere_scale, 0.4